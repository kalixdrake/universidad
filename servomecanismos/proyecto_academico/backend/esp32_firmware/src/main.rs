//! Firmware ESP32 — Robot 2R
//!
//! Control PID de 2 motores DC con encoders magnéticos (16 PPR, cuadratura ×4)
//! y puentes H BTS7960 (base) + L298N (codo).
//!
//! ## Parámetros de motores (validados contra modelo_motor_dc.m)
//!
//! ### Motor Base — 12 V, N = 218:1
//!   R = 1.200 Ω,  Ke = 0.0142 V·s/rad,  Kt = 0.0142 N·m/A
//!   ω_nl (motor) = 844.7 rad/s  (8066 RPM)
//!   ω_nl (salida) = 3.87 rad/s   (37 RPM)
//!   Encoder efectivo = 4 × 16 × 218 = 13 952 cuentas/rev @ salida
//!   Resolución angular @ salida = 0.0258° por cuenta
//!
//! ### Motor Codo — 24 V, N = 270:1
//!   R = 7.742 Ω,  Ke = 0.0229 V·s/rad,  Kt = 0.0229 N·m/A
//!   ω_nl (motor) = 1046.2 rad/s  (9990 RPM)
//!   ω_nl (salida) = 3.87 rad/s   (37 RPM)
//!   Encoder efectivo = 4 × 16 × 270 = 17 280 cuentas/rev @ salida
//!   Resolución angular @ salida = 0.0208° por cuenta
//!
//! ## Frecuencias
//!   PWM:     25 kHz  (MCPWM, autónomo de CPU)
//!   PID:      1 kHz  (~5 cuentas encoder por iteración)
//!   Telemetría: 100 Hz  (~53 cuentas, ~10 iteraciones PID por frame)
//!
//! ## Pines
//!   Base (BTS7960):  RPWM=GPIO25, LPWM=GPIO26
//!   Codo (L298N):    IN1=GPIO14, IN2=GPIO27, ENA=GPIO12
//!   Encoder Base:    A=GPIO34, B=GPIO35  (PCNT unit 0)
//!   Encoder Codo:    A=GPIO36, B=GPIO39  (PCNT unit 1)
//!   UART0 (USB):     TX=GPIO1, RX=GPIO3

#![no_std]
#![no_main]

use esp_backtrace as _;
use esp_hal::{
    clock::ClockControl,
    peripherals::Peripherals,
    prelude::*,
    Delay,
    Uart,
};
use esp_println::println;

// ═══════════════════════════════════════════════════════════════
//  PARÁMETROS DE MOTORES (del modelo_motor_dc.m)
// ═══════════════════════════════════════════════════════════════

/// Motor Base — JGB37-3530, 12 V, 30 kgf·cm nominal
mod motor_base {
    pub const V_NOM: f32 = 12.0;           // [V]
    pub const R: f32 = 1.200;              // [Ω]
    pub const KE: f32 = 0.0142;            // [V·s/rad] back-EMF
    pub const KT: f32 = 0.0142;            // [N·m/A] constante de torque (=KE)
    pub const N_REDUCTION: f32 = 218.0;    // Relación caja planetaria
    pub const J_ROTOR: f32 = 3.5e-6;       // [kg·m²] inercia del rotor
    pub const ENCODER_PPR: u32 = 16;       // Pulsos por revolución (eje motor)
    // Efectivo a la salida: 4 × PPR × N = 13 952 cuentas/rev
    pub const COUNTS_PER_REV_OUT: f32 = 4.0 * ENCODER_PPR as f32 * N_REDUCTION;
    pub const RAD_PER_COUNT: f32 = core::f32::consts::TAU / COUNTS_PER_REV_OUT;
    pub const OMEGA_MAX: f32 = 844.7;       // ω sin carga @ 12V [rad/s]
}

/// Motor Codo — JGB37-3530, 24 V, 19 kgf·cm nominal
mod motor_codo {
    pub const V_NOM: f32 = 24.0;
    pub const R: f32 = 7.742;
    pub const KE: f32 = 0.0229;
    pub const KT: f32 = 0.0229;
    pub const N_REDUCTION: f32 = 270.0;
    pub const J_ROTOR: f32 = 4.0e-6;
    pub const ENCODER_PPR: u32 = 16;
    pub const COUNTS_PER_REV_OUT: f32 = 4.0 * ENCODER_PPR as f32 * N_REDUCTION;
    pub const RAD_PER_COUNT: f32 = core::f32::consts::TAU / COUNTS_PER_REV_OUT;
    pub const OMEGA_MAX: f32 = 1046.2;
}

// ═══════════════════════════════════════════════════════════════
//  FRECUENCIAS DEL SISTEMA
// ═══════════════════════════════════════════════════════════════
const PWM_FREQ_HZ: u32 = 25_000;
const CONTROL_FREQ_HZ: u32 = 1_000;
const CONTROL_PERIOD_US: u64 = 1_000_000 / CONTROL_FREQ_HZ as u64;
const TELEMETRY_FREQ_HZ: u32 = 100;
const TELEMETRY_PERIOD_US: u64 = 1_000_000 / TELEMETRY_FREQ_HZ as u64;

// ═══════════════════════════════════════════════════════════════
//  GANANCIAS PID  (sintonizables en runtime vía comandos)
// ═══════════════════════════════════════════════════════════════

// Base: τ_m ≈ 21 ms → polos en ~50 rad/s → Kp ≈ 3.0, Ki ≈ 0.15
const PID_BASE_KP: f32 = 3.0;
const PID_BASE_KI: f32 = 0.15;
const PID_BASE_KD: f32 = 0.02;
const PID_BASE_I_MAX: f32 = 0.8;

// Codo: τ_m ≈ 59 ms → polos en ~17 rad/s → Kp ≈ 2.5, Ki ≈ 0.10
const PID_CODO_KP: f32 = 2.5;
const PID_CODO_KI: f32 = 0.10;
const PID_CODO_KD: f32 = 0.015;
const PID_CODO_I_MAX: f32 = 0.8;

// ═══════════════════════════════════════════════════════════════
//  ESTRUCTURAS DE DATOS
// ═══════════════════════════════════════════════════════════════

const SYNC_WORD: [u8; 2] = [0xAA, 0x55];
const MAX_WAYPOINTS: usize = 4096;

#[derive(Clone, Copy)]
struct Waypoint {
    theta1: f32,   // [rad] posición deseada a la salida
    theta2: f32,
    omega1: f32,   // [rad/s] velocidad deseada a la salida
    omega2: f32,
    dt: f32,       // [s] paso de tiempo
}

#[derive(Clone, Copy, PartialEq)]
enum RobotState {
    Idle,
    Receiving,
    Running,
    Error(u8),
    Completed,
}

/// Controlador PID con anti-windup + feedforward de back-EMF.
struct MotorController {
    kp: f32,
    ki: f32,
    kd: f32,
    integral: f32,
    i_max: f32,
    prev_error: f32,
    ke: f32,            // [V·s/rad] back-EMF
    v_nom: f32,         // [V] voltaje nominal
    output_min: f32,
    output_max: f32,
    encoder_counts: i32,
    encoder_prev: i32,
    vel_samples: [f32; 4],
    vel_idx: usize,
}

impl MotorController {
    fn new(kp: f32, ki: f32, kd: f32, i_max: f32, ke: f32, v_nom: f32) -> Self {
        Self {
            kp, ki, kd,
            integral: 0.0, i_max,
            prev_error: 0.0,
            ke, v_nom,
            output_min: -1.0, output_max: 1.0,
            encoder_counts: 0, encoder_prev: 0,
            vel_samples: [0.0; 4], vel_idx: 0,
        }
    }

    fn counts_to_rad(counts: i32, rpc: f32) -> f32 {
        counts as f32 * rpc
    }

    /// Estima velocidad [rad/s] con media móvil de 4 muestras.
    fn update_velocity(&mut self, new_counts: i32, dt: f32, rpc: f32) -> f32 {
        let delta = new_counts.wrapping_sub(self.encoder_prev) as f32;
        self.encoder_prev = new_counts;
        self.encoder_counts = new_counts;
        let inst = delta * rpc / dt;
        self.vel_samples[self.vel_idx] = inst;
        self.vel_idx = (self.vel_idx + 1) % 4;
        self.vel_samples.iter().sum::<f32>() / 4.0
    }

    /// PID + feedforward.
    /// sp_pos, meas_pos: [rad]  |  sp_vel, meas_vel: [rad/s]  |  dt: [s]
    fn compute(&mut self, sp_pos: f32, sp_vel: f32, meas_pos: f32, meas_vel: f32, dt: f32) -> f32 {
        let error = sp_pos - meas_pos;
        self.integral = (self.integral + error * dt).clamp(-self.i_max, self.i_max);
        let derivative = sp_vel - meas_vel; // derivada suave (evita amplificar ruido de encoder)
        let pid = self.kp * error + self.ki * self.integral + self.kd * derivative;
        let ff = self.ke * sp_vel / self.v_nom; // feedforward back-EMF
        self.prev_error = error;
        (pid + ff).clamp(self.output_min, self.output_max)
    }

    fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.vel_samples = [0.0; 4];
        self.vel_idx = 0;
    }
}

// ═══════════════════════════════════════════════════════════════
//  BUFFERS GLOBALES
// ═══════════════════════════════════════════════════════════════

static mut TRAJECTORY: [Waypoint; MAX_WAYPOINTS] = [Waypoint {
    theta1: 0.0, theta2: 0.0, omega1: 0.0, omega2: 0.0, dt: 0.01,
}; MAX_WAYPOINTS];
static mut TRAJECTORY_LEN: usize = 0;
static mut WAYPOINT_IDX: usize = 0;
static mut STATE: RobotState = RobotState::Idle;
static mut CYCLE_COUNT: u32 = 1;
static mut CYCLE_DONE: u32 = 0;

// ═══════════════════════════════════════════════════════════════
//  PROTOCOLO (compatible con backend/communication/)
// ═══════════════════════════════════════════════════════════════

mod protocol {
    use serde::{Deserialize, Serialize};
    #[derive(Debug, Serialize, Deserialize)]
    pub enum Frame {
        NetworkConfig { ssid: heapless::String<32>, password: heapless::String<64>, static_ip: heapless::String<16> },
        TrajectoryMeta { dt: f32, cycles: u32, total_points: u32 },
        TrajectoryChunk {
            chunk_index: u32, total_chunks: u32, num_points: u32,
            theta1: heapless::Vec<f32, 256>, theta2: heapless::Vec<f32, 256>,
            omega1: heapless::Vec<f32, 256>, omega2: heapless::Vec<f32, 256>,
        },
        StartExecution,
        EmergencyStop,
        Ack { success: bool, message: heapless::String<64> },
        Telemetry { timestamp_us: u64, theta1: f32, theta2: f32, omega1: f32, omega2: f32, waypoint_index: u32, state: u8, error_code: u8 },
        Heartbeat { seq: u32 },
    }
}

// ═══════════════════════════════════════════════════════════════
//  PUNTO DE ENTRADA
// ═══════════════════════════════════════════════════════════════

#[entry]
fn main() -> ! {
    let peripherals = Peripherals::take();
    let system = peripherals.SYSTEM.split();
    let clocks = ClockControl::max(system.clock_control).freeze();
    let mut delay = Delay::new(&clocks);
    let mut uart0 = Uart::new(peripherals.UART0, &clocks);

    println!("╔══════════════════════════════════════╗");
    println!("║   ESP32 Robot 2R Firmware v1.0      ║");
    println!("╠══════════════════════════════════════╣");
    println!("║ Base: 12V N=218 16PPR  BTS7960      ║");
    println!("║ Codo: 24V N=270 16PPR  L298N        ║");
    println!("║ PID: {} Hz | Telem: {} Hz         ║", CONTROL_FREQ_HZ, TELEMETRY_FREQ_HZ);
    println!("║ PWM: {} kHz | MCPWM autónomo     ║", PWM_FREQ_HZ / 1000);
    println!("╚══════════════════════════════════════╝");

    // ── Controladores ──────────────────────────────
    let mut ctrl_base = MotorController::new(
        PID_BASE_KP, PID_BASE_KI, PID_BASE_KD,
        PID_BASE_I_MAX, motor_base::KE, motor_base::V_NOM,
    );
    let mut ctrl_codo = MotorController::new(
        PID_CODO_KP, PID_CODO_KI, PID_CODO_KD,
        PID_CODO_I_MAX, motor_codo::KE, motor_codo::V_NOM,
    );

    // ── Variables de tiempo ────────────────────────
    let dt = CONTROL_PERIOD_US as f32 / 1_000_000.0;
    let mut micros: u64 = 0;
    let mut last_telem_us: u64 = 0;

    println!("Firmware listo. Esperando trayectoria...\n");

    // ═══════════════════════════════════════════════════
    //  BUCLE PRINCIPAL  (1 kHz)
    // ═══════════════════════════════════════════════════
    loop {
        micros = micros.wrapping_add(CONTROL_PERIOD_US);

        // ── 1. Leer encoders (PCNT por hardware) ─────
        // TODO: pcnt0.get_value(), pcnt1.get_value()
        let raw_base: i32 = 0;
        let raw_codo: i32 = 0;

        let pos_base = MotorController::counts_to_rad(raw_base, motor_base::RAD_PER_COUNT);
        let pos_codo = MotorController::counts_to_rad(raw_codo, motor_codo::RAD_PER_COUNT);
        let vel_base = ctrl_base.update_velocity(raw_base, dt, motor_base::RAD_PER_COUNT);
        let vel_codo = ctrl_codo.update_velocity(raw_codo, dt, motor_codo::RAD_PER_COUNT);

        // ── 2. Control PID ───────────────────────────
        let (d1, d2) = unsafe {
            match STATE {
                RobotState::Running => {
                    if WAYPOINT_IDX < TRAJECTORY_LEN {
                        let wp = &TRAJECTORY[WAYPOINT_IDX];
                        let duty1 = ctrl_base.compute(wp.theta1, wp.omega1, pos_base, vel_base, wp.dt);
                        let duty2 = ctrl_codo.compute(wp.theta2, wp.omega2, pos_codo, vel_codo, wp.dt);
                        WAYPOINT_IDX += 1;
                        (duty1, duty2)
                    } else {
                        CYCLE_DONE += 1;
                        if CYCLE_DONE >= CYCLE_COUNT {
                            STATE = RobotState::Completed;
                            ctrl_base.reset();
                            ctrl_codo.reset();
                        } else {
                            WAYPOINT_IDX = 0;
                        }
                        (0.0, 0.0)
                    }
                }
                _ => (0.0, 0.0),
            }
        };

        // ── 3. Aplicar PWM (MCPWM autónomo) ─────────
        set_motor_base(d1);
        set_motor_codo(d2);

        // ── 4. Telemetría a 100 Hz ───────────────────
        if micros.wrapping_sub(last_telem_us) >= TELEMETRY_PERIOD_US {
            last_telem_us = micros;
            send_telemetry(&mut uart0, micros, pos_base, pos_codo, vel_base, vel_codo);
        }

        // ── 5. Esperar fin de período ────────────────
        delay.delay_micros((CONTROL_PERIOD_US - 80) as u32);
    }
}

// ═══════════════════════════════════════════════════════════════
//  HARDWARE: PWM MOTORES
// ═══════════════════════════════════════════════════════════════

/// BTS7960 (Base): RPWM=GPIO25, LPWM=GPIO26.
/// duty > 0 → RPWM=|duty|, LPWM=0.  duty < 0 → RPWM=0, LPWM=|duty|.
fn set_motor_base(duty: f32) {
    // TODO: MCPWM0 chA=GPIO25, chB=GPIO26 @ 25 kHz
    let _ = duty;
}

/// L298N (Codo): IN1=GPIO14, IN2=GPIO27, ENA=GPIO12.
/// duty > 0 → IN1=H, IN2=L, ENA=|duty|.  duty < 0 → IN1=L, IN2=H, ENA=|duty|.
fn set_motor_codo(duty: f32) {
    // TODO: GPIO14/27 (dir), MCPWM1 chA=GPIO12 (ENA) @ 25 kHz
    let _ = duty;
}

// ═══════════════════════════════════════════════════════════════
//  TELEMETRÍA
// ═══════════════════════════════════════════════════════════════

fn send_telemetry(
    _uart: &mut Uart<'_, esp_hal::peripherals::UART0>,
    ts: u64, th1: f32, th2: f32, om1: f32, om2: f32,
) {
    let (st, wp) = unsafe {
        let code: u8 = match STATE {
            RobotState::Idle | RobotState::Receiving => 0,
            RobotState::Running => 1,
            RobotState::Error(_) => 2,
            RobotState::Completed => 3,
        };
        (code, WAYPOINT_IDX as u32)
    };

    let mut payload = [0u8; 30];
    payload[0..8].copy_from_slice(&ts.to_le_bytes());
    payload[8..12].copy_from_slice(&th1.to_le_bytes());
    payload[12..16].copy_from_slice(&th2.to_le_bytes());
    payload[16..20].copy_from_slice(&om1.to_le_bytes());
    payload[20..24].copy_from_slice(&om2.to_le_bytes());
    payload[24..28].copy_from_slice(&wp.to_le_bytes());
    payload[28] = st;
    payload[29] = 0;

    let mut frame = [0u8; 39];
    frame[0..2].copy_from_slice(&SYNC_WORD);
    frame[2..4].copy_from_slice(&(30u16).to_le_bytes());
    frame[4] = 0x30;
    frame[5..35].copy_from_slice(&payload);
    let crc = crc32_soft(&frame[..35]);
    frame[35..39].copy_from_slice(&crc.to_le_bytes());

    // TODO: uart.write_bytes(&frame)
    let _ = (ts, th1, th2, om1, om2);
}

/// CRC32 software (compatible con crc32fast de rust_comm).
fn crc32_soft(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFF_FFFF;
    for &b in data {
        crc ^= b as u32;
        for _ in 0..8 {
            crc = if crc & 1 != 0 { (crc >> 1) ^ 0xEDB8_8320 } else { crc >> 1 };
        }
    }
    !crc
}
