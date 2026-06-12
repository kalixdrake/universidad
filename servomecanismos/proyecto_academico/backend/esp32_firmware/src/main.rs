//! Firmware para ESP32 del robot 2R.
//!
//! ## Funcionalidad
//! - Recibe configuración de red vía USB/serial
//! - Almacena configuración WiFi en NVS
//! - Recibe trayectorias (θ1, θ2, ω1, ω2) vía TCP o serial
//! - Ejecuta control PID en tiempo real sobre dos motores DC
//! - Envía telemetría (posición, velocidad) a 100 Hz
//!
//! ## Arquitectura
//! ```
//! ┌─────────────────────────────────────────────┐
//! │  main loop (1 kHz)                          │
//! │  ├── Lee encoders (QEI)                     │
//! │  ├── Calcula PID                             │
//! │  ├── Actualiza PWM                           │
//! │  └── Cada 10 ms: envía telemetría            │
//! │                                              │
//! │  Tarea WiFi (background)                     │
//! │  ├── Acepta conexiones TCP                   │
//! │  ├── Recibe chunks de trayectoria            │
//! │  └── Responde a UDP broadcast discovery       │
//! └─────────────────────────────────────────────┘
//! ```

#![no_std]
#![no_main]

use esp_backtrace as _;
use esp_hal::{
    clock::ClockControl,
    gpio::IO,
    peripherals::Peripherals,
    prelude::*,
    timer::TimerGroup,
    Delay,
    Rtc,
    Uart,
};
use esp_println::println;

// ── Constantes del robot ──────────────────────────────
const CONTROL_FREQ_HZ: u32 = 1000;
const TELEMETRY_FREQ_HZ: u32 = 100;
const PID_KP: f32 = 2.5;
const PID_KI: f32 = 0.1;
const PID_KD: f32 = 0.05;
const ENCODER_CPR: u32 = 360; // pulsos por revolución
const SYNC_WORD: [u8; 2] = [0xAA, 0x55];

// ── Estructuras de datos ──────────────────────────────

/// Configuración de red almacenada en NVS.
#[derive(serde::Serialize, serde::Deserialize, Clone)]
struct StoredConfig {
    ssid: heapless::String<32>,
    password: heapless::String<64>,
    static_ip: heapless::String<16>,
}

/// Punto de trayectoria.
#[derive(Clone, Copy)]
struct Waypoint {
    theta1: f32,
    theta2: f32,
    omega1: f32,
    omega2: f32,
    dt: f32,
}

/// Estado del robot.
#[derive(Clone, Copy, PartialEq)]
enum RobotState {
    Idle,
    Receiving,
    Running,
    Error(u8),
    Completed,
}

/// Controlador PID simple.
struct PidController {
    kp: f32,
    ki: f32,
    kd: f32,
    integral: f32,
    prev_error: f32,
    setpoint: f32,
    output_min: f32,
    output_max: f32,
}

impl PidController {
    fn new(kp: f32, ki: f32, kd: f32, out_min: f32, out_max: f32) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            setpoint: 0.0,
            output_min: out_min,
            output_max: out_max,
        }
    }

    fn set_setpoint(&mut self, sp: f32) {
        self.setpoint = sp;
    }

    fn update(&mut self, measurement: f32, dt: f32) -> f32 {
        let error = self.setpoint - measurement;
        self.integral += error * dt;
        let derivative = (error - self.prev_error) / dt.max(1e-6);
        self.prev_error = error;

        let output = self.kp * error + self.ki * self.integral + self.kd * derivative;
        output.clamp(self.output_min, self.output_max)
    }

    fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
}

// ── Buffers globales ──────────────────────────────────
static mut TRAJECTORY_BUFFER: [Waypoint; 4096] = [Waypoint {
    theta1: 0.0,
    theta2: 0.0,
    omega1: 0.0,
    omega2: 0.0,
    dt: 0.01,
}; 4096];
static mut TRAJECTORY_LEN: usize = 0;
static mut CURRENT_WAYPOINT: usize = 0;
static mut ROBOT_STATE: RobotState = RobotState::Idle;
static mut CYCLE_COUNT: u32 = 0;
static mut CYCLES_COMPLETED: u32 = 0;

// ── Protocolo (mismo que en PC) ───────────────────────
mod protocol {
    use serde::{Deserialize, Serialize};

    #[derive(Debug, Serialize, Deserialize)]
    pub enum Frame {
        NetworkConfig {
            ssid: heapless::String<32>,
            password: heapless::String<64>,
            static_ip: heapless::String<16>,
        },
        TrajectoryChunk {
            chunk_index: u32,
            total_chunks: u32,
            theta1: heapless::Vec<f32, 256>,
            theta2: heapless::Vec<f32, 256>,
            omega1: heapless::Vec<f32, 256>,
            omega2: heapless::Vec<f32, 256>,
            dt: f32,
            cycles: u32,
        },
        StartExecution,
        EmergencyStop,
        Ack {
            success: bool,
            message: heapless::String<64>,
        },
        Telemetry {
            timestamp_us: u64,
            theta1: f32,
            theta2: f32,
            omega1: f32,
            omega2: f32,
            waypoint_index: u32,
            state: u8,
            error_code: u8,
        },
        Heartbeat {
            seq: u32,
        },
    }
}

// ── Punto de entrada ──────────────────────────────────
#[entry]
fn main() -> ! {
    let peripherals = Peripherals::take();
    let system = peripherals.SYSTEM.split();
    let clocks = ClockControl::max(system.clock_control).freeze();
    let mut delay = Delay::new(&clocks);

    let io = IO::new(peripherals.GPIO, peripherals.IO_MUX);
    let mut uart0 = Uart::new(peripherals.UART0, &clocks);

    println!("=== ESP32 Robot 2R Firmware v0.1 ===");
    println!("Esperando configuración por USB...");

    // TODO: Inicializar WiFi, PWM, encoders
    // TODO: Cargar StoredConfig desde NVS
    // TODO: Iniciar tarea de descubrimiento UDP

    let mut pid_motor1 = PidController::new(PID_KP, PID_KI, PID_KD, -1.0, 1.0);
    let mut pid_motor2 = PidController::new(PID_KP, PID_KI, PID_KD, -1.0, 1.0);

    let mut last_telemetry = 0u64;
    let mut micros: u64 = 0;

    println!("Firmware listo. Esperando datos...");

    loop {
        micros += 1000; // 1 kHz -> 1000 µs por iteración

        // ── Lectura de comandos por serial ────────────
        // (implementación del protocolo de tramas)

        // ── Control PID ───────────────────────────────
        unsafe {
            if ROBOT_STATE == RobotState::Running {
                if CURRENT_WAYPOINT < TRAJECTORY_LEN {
                    let wp = &TRAJECTORY_BUFFER[CURRENT_WAYPOINT];

                    pid_motor1.set_setpoint(wp.theta1);
                    pid_motor2.set_setpoint(wp.theta2);

                    // TODO: Leer encoders reales
                    let enc1: f32 = 0.0;
                    let enc2: f32 = 0.0;

                    let pwm1 = pid_motor1.update(enc1, wp.dt);
                    let pwm2 = pid_motor2.update(enc2, wp.dt);

                    // TODO: Aplicar PWM a motores

                    // Avanzar al siguiente waypoint
                    CURRENT_WAYPOINT += 1;
                } else {
                    CYCLES_COMPLETED += 1;
                    if CYCLES_COMPLETED >= CYCLE_COUNT {
                        ROBOT_STATE = RobotState::Completed;
                    } else {
                        CURRENT_WAYPOINT = 0; // reiniciar ciclo
                    }
                }
            }
        }

        // ── Envío de telemetría ───────────────────────
        if micros - last_telemetry >= 1_000_000 / TELEMETRY_FREQ_HZ as u64 {
            last_telemetry = micros;
            // TODO: Enviar Frame::Telemetry por UART/WiFi
            unsafe {
                let state_code: u8 = match ROBOT_STATE {
                    RobotState::Idle => 0,
                    RobotState::Receiving => 0,
                    RobotState::Running => 1,
                    RobotState::Error(code) => 2 + code,
                    RobotState::Completed => 3,
                };
                // Serializar y enviar
            }
        }

        delay.delay_micros(1000 - 50); // compensar overhead ~50 µs
    }
}
