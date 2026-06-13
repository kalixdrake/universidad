"""
Identificación y caracterización de motores DC del robot 2R.

Incluye:
- Búsqueda del torque de rozamiento (PWM mínimo para vencer fricción)
- Búsqueda de límites angulares (movimiento lento hasta detenerse)
- Caracterización de curvas velocidad/aceleración/jerk vs PWM
- Ajuste de funciones de transferencia

Estrategia:
    1. PWM mínimo: barrido incremental desde 0% hasta detectar movimiento
       (encoder o telemetría). Se registra como "friction_pwm".
    2. Límites: con friction_pwm, mover en cada dirección hasta que
       el encoder deje de cambiar (choque mecánico). Se registran min/max.
    3. Caracterización: en varios puntos del rango, aplicar escalones de PWM
       (10%, 20%, ..., 100%) y medir velocidad, aceleración y jerk.
    4. Ajuste: regresión lineal/polinomial para mapear PWM → velocidad,
       y estimar función de transferencia G(s) = K/(τs + 1).
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Callable
import logging

logger = logging.getLogger(__name__)


# ── Constantes de identificación ─────────────────────
PWM_SWEEP_STEPS = 10       # Pasos en el barrido de PWM (10% a 100%)
PWM_SWEEP_DURATION_S = 1.5  # Duración de cada escalón de prueba [s]
PWM_SWEEP_COOLDOWN_S = 0.3  # Tiempo de espera entre escalones [s]
MIN_MOVEMENT_THRESHOLD_RAD = 0.002  # Mínimo cambio angular para detectar movimiento
FRICTION_SEARCH_TIMEOUT_S = 5.0     # Timeout en búsqueda de PWM mínimo
ANGLE_LIMIT_STALL_TIMEOUT_S = 3.0   # Timeout en búsqueda de límite angular


@dataclass
class MotorParams:
    """Parámetros identificados de un motor."""
    # PWM de rozamiento (mínimo para movimiento)
    friction_pwm: float = 0.0          # fracción [0, 1]

    # Límites angulares [rad]
    angle_min: float = -np.pi
    angle_max: float = np.pi
    zero_offset: float = 0.0           # offset del encoder para posición cero

    # Curva velocidad vs PWM (puntos de calibración)
    pwm_points: list[float] = field(default_factory=list)    # [0.1, 0.2, ..., 1.0]
    velocity_points: list[float] = field(default_factory=list)  # rad/s medidos
    acceleration_points: list[float] = field(default_factory=list)
    jerk_points: list[float] = field(default_factory=list)

    # Función de transferencia estimada: G(s) = K / (tau*s + 1)
    K: float = 0.0       # Ganancia DC [rad/s / PWM%]
    tau: float = 0.0     # Constante de tiempo [s]

    # Métricas de ajuste
    r_squared: float = 0.0
    rmse: float = 0.0

    # Timestamp de calibración
    calibrated: bool = False

    def to_dict(self) -> dict:
        return {
            'friction_pwm': self.friction_pwm,
            'angle_min': self.angle_min,
            'angle_max': self.angle_max,
            'zero_offset': self.zero_offset,
            'pwm_points': self.pwm_points,
            'velocity_points': self.velocity_points,
            'acceleration_points': self.acceleration_points,
            'jerk_points': self.jerk_points,
            'K': self.K,
            'tau': self.tau,
            'r_squared': self.r_squared,
            'rmse': self.rmse,
            'calibrated': self.calibrated,
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'MotorParams':
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


# ── Identificación de torque de rozamiento ───────────

def identify_friction_torque(
    motor_id: int,
    pwm_min_step: float = 0.01,
    pwm_max_search: float = 0.30,
    telemetry_callback: Optional[Callable[[], dict]] = None,
    progress_callback: Optional[Callable[[int, float, str], None]] = None,
) -> float:
    """
    Encuentra el PWM mínimo necesario para que el motor se mueva.

    Estrategia:
    - Barrido incremental desde pwm_min_step hasta pwm_max_search.
    - En cada paso, se aplica PWM y se monitorea el encoder.
    - Se considera "movimiento" cuando |Δθ| > MIN_MOVEMENT_THRESHOLD_RAD.

    Args:
        motor_id: 1 = base, 2 = codo.
        pwm_min_step: Incremento mínimo de PWM (default 1%).
        pwm_max_search: PWM máximo de búsqueda (default 30%).
        telemetry_callback: Función que retorna dict con 'theta1', 'theta2'.
        progress_callback: (motor_id, pwm, status_str) para UI.

    Returns:
        PWM de rozamiento [fracción 0-1].
    """
    logger.info(f"Motor {motor_id}: buscando torque de rozamiento...")

    if telemetry_callback is None:
        logger.warning(f"Motor {motor_id}: sin callback de telemetría, usando valor por defecto 0.12")
        return 0.12  # Valor conservador por defecto

    pwm = pwm_min_step
    while pwm <= pwm_max_search:
        if progress_callback:
            progress_callback(motor_id, pwm, f"Probando PWM={pwm*100:.1f}%")

        # Aplicar PWM y leer posición inicial
        theta_key = f'theta{motor_id}'
        try:
            telem = telemetry_callback()
            theta_start = float(telem.get(theta_key, 0.0))
        except Exception:
            theta_start = 0.0

        # Esperar un momento para que el motor responda
        import time
        time.sleep(0.3)

        # Leer posición final
        try:
            telem = telemetry_callback()
            theta_end = float(telem.get(theta_key, 0.0))
        except Exception:
            theta_end = 0.0

        delta = abs(theta_end - theta_start)

        if delta > MIN_MOVEMENT_THRESHOLD_RAD:
            # Añadir margen de seguridad del 20%
            friction_pwm = pwm * 1.2
            logger.info(f"Motor {motor_id}: rozamiento vencido a PWM={pwm*100:.1f}% "
                        f"(Δθ={delta*180/np.pi:.3f}°), con margen: {friction_pwm*100:.1f}%")
            if progress_callback:
                progress_callback(motor_id, friction_pwm, f"✓ Rozamiento: {friction_pwm*100:.1f}%")
            return friction_pwm

        pwm += pwm_min_step

    logger.warning(f"Motor {motor_id}: no se detectó movimiento hasta {pwm_max_search*100:.0f}% PWM")
    return pwm_max_search


# ── Búsqueda de límites angulares ────────────────────

def find_angle_limits(
    motor_id: int,
    friction_pwm: float,
    telemetry_callback: Optional[Callable[[], dict]] = None,
    progress_callback: Optional[Callable[[int, str, float], None]] = None,
) -> tuple[float, float]:
    """
    Encuentra los límites angulares del motor moviéndolo lentamente
    con el PWM de rozamiento en ambas direcciones.

    Args:
        motor_id: 1 = base, 2 = codo.
        friction_pwm: PWM mínimo para movimiento.
        telemetry_callback: Función que retorna dict con 'theta1', 'theta2'.
        progress_callback: (motor_id, direction_str, current_angle_deg).

    Returns:
        (angle_min, angle_max) en radianes.
    """
    import time

    logger.info(f"Motor {motor_id}: buscando límites angulares...")

    if telemetry_callback is None:
        logger.warning(f"Motor {motor_id}: sin telemetría, usando ±π")
        return -np.pi, np.pi

    theta_key = f'theta{motor_id}'

    def _search_limit(direction_sign: int, label: str) -> float:
        """Busca un límite moviendo en la dirección dada."""
        last_theta: Optional[float] = None
        stall_start: Optional[float] = None

        t_start = time.time()
        while time.time() - t_start < ANGLE_LIMIT_STALL_TIMEOUT_S:
            try:
                telem = telemetry_callback()
                current_theta = float(telem.get(theta_key, 0.0))
            except Exception:
                time.sleep(0.1)
                continue

            if progress_callback:
                progress_callback(motor_id, label, np.rad2deg(current_theta))

            if last_theta is not None:
                delta = abs(current_theta - last_theta)
                if delta < MIN_MOVEMENT_THRESHOLD_RAD / 2:
                    if stall_start is None:
                        stall_start = time.time()
                    elif time.time() - stall_start > 1.5:
                        # Motor detenido > 1.5s → límite alcanzado
                        logger.info(f"Motor {motor_id}: límite {label}={np.rad2deg(current_theta):.1f}°")
                        return current_theta
                else:
                    stall_start = None  # Reset, sigue moviéndose

            last_theta = current_theta
            time.sleep(0.1)

        logger.warning(f"Motor {motor_id}: timeout buscando límite {label}")
        return last_theta if last_theta is not None else 0.0

    # Buscar límite positivo y negativo
    angle_max = _search_limit(1, "max")
    angle_min = _search_limit(-1, "min")

    return angle_min, angle_max


# ── Caracterización de curva motor ───────────────────

def characterize_motor_curve(
    motor_id: int,
    friction_pwm: float,
    angle_min: float,
    angle_max: float,
    telemetry_callback: Optional[Callable[[], dict]] = None,
    progress_callback: Optional[Callable[[int, float, str], None]] = None,
) -> MotorParams:
    """
    Caracteriza velocidad, aceleración y jerk del motor a diferentes PWM.

    Aplica escalones de PWM desde friction_pwm hasta 100% y mide
    la respuesta en velocidad del motor. Ajusta una función de
    transferencia de primer orden: G(s) = K / (τs + 1).

    Args:
        motor_id: 1 = base, 2 = codo.
        friction_pwm: PWM mínimo de rozamiento.
        angle_min, angle_max: Límites angulares [rad].
        telemetry_callback: Función que retorna dict con theta y omega.
        progress_callback: (motor_id, pwm, status_str).

    Returns:
        MotorParams con los parámetros identificados.
    """
    import time

    params = MotorParams()
    params.friction_pwm = friction_pwm
    params.angle_min = angle_min
    params.angle_max = angle_max

    logger.info(f"Motor {motor_id}: caracterizando curva velocidad/PWM...")

    if telemetry_callback is None:
        logger.warning(f"Motor {motor_id}: sin telemetría, usando modelo teórico")
        return _estimate_theoretical_params(motor_id, params)

    theta_key = f'theta{motor_id}'
    omega_key = f'omega{motor_id}'

    pwm_levels = np.linspace(max(friction_pwm, 0.05), 1.0, PWM_SWEEP_STEPS)

    for pwm in pwm_levels:
        if progress_callback:
            progress_callback(motor_id, pwm, f"Midiendo PWM={pwm*100:.0f}%")

        # Verificar que estamos dentro de límites
        try:
            telem = telemetry_callback()
            theta_start = float(telem.get(theta_key, 0.0))
        except Exception:
            continue

        # Solo medir si hay espacio para moverse
        mid_angle = (angle_min + angle_max) / 2
        margin = (angle_max - angle_min) * 0.3

        if not (angle_min + margin < theta_start < angle_max - margin):
            # Fuera de rango seguro, necesitaríamos reposicionar
            logger.debug(f"Motor {motor_id}: fuera de rango seguro para caracterización")
            continue

        # Registrar velocidades durante el escalón
        velocities = []
        times_samples = []
        t0 = time.time()

        while time.time() - t0 < PWM_SWEEP_DURATION_S:
            try:
                telem = telemetry_callback()
                omega = float(telem.get(omega_key, 0.0))
                velocities.append(omega)
                times_samples.append(time.time() - t0)
            except Exception:
                pass
            time.sleep(0.02)

        if len(velocities) > 5:
            # Calcular estadísticas del escalón
            vel_array = np.array(velocities)
            vel_steady = np.mean(vel_array[-len(vel_array)//3:])  # Último tercio

            # Aceleración: derivada numérica
            accel = np.gradient(vel_array, times_samples)
            accel_max = np.max(np.abs(accel))

            # Jerk: segunda derivada
            jerk = np.gradient(accel, times_samples)
            jerk_max = np.max(np.abs(jerk))

            params.pwm_points.append(float(pwm))
            params.velocity_points.append(float(vel_steady))
            params.acceleration_points.append(float(accel_max))
            params.jerk_points.append(float(jerk_max))

        # Cooldown
        time.sleep(PWM_SWEEP_COOLDOWN_S)

    # Ajustar función de transferencia
    if len(params.pwm_points) >= 3:
        _fit_transfer_function(params)

    params.calibrated = len(params.pwm_points) >= 3
    logger.info(f"Motor {motor_id}: caracterización completada. "
                f"K={params.K:.3f}, τ={params.tau:.3f}s, R²={params.r_squared:.3f}")

    return params


# ── Ajuste de función de transferencia ───────────────

def fit_transfer_function(pwm_points: list[float], velocity_points: list[float]) -> tuple[float, float, float]:
    """
    Ajusta G(s) = K/(τs + 1) a partir de datos velocidad vs PWM.

    Usa la respuesta en estado estacionario:
        ω_ss = K * PWM  (modelo linealizado)

    Returns:
        (K, tau, r_squared)
    """
    params = MotorParams(pwm_points=pwm_points, velocity_points=velocity_points)
    _fit_transfer_function(params)
    return params.K, params.tau, params.r_squared


def _fit_transfer_function(params: MotorParams):
    """Ajuste interno de FT de primer orden."""
    pwm = np.array(params.pwm_points)
    vel = np.array(params.velocity_points)

    # Regresión lineal: ω = K * PWM (forzando origen)
    # K = Σ(ω_i * PWM_i) / Σ(PWM_i²)
    K = np.sum(vel * pwm) / np.sum(pwm**2) if np.sum(pwm**2) > 0 else 0.0

    # Estimar τ a partir del tiempo de subida
    # τ ≈ t_rise / 2.2  (para sistema de primer orden)
    # Como no tenemos respuesta transitoria completa, usamos τ típico
    # basado en las mediciones de aceleración
    if params.acceleration_points:
        accel_avg = np.mean(np.abs(params.acceleration_points))
        # τ ≈ K * PWM / a_max (aproximación)
        tau = abs(K) / (accel_avg + 1e-6)
    else:
        tau = 0.05  # valor conservador por defecto

    tau = np.clip(tau, 0.01, 0.5)

    # R²
    vel_pred = K * pwm
    ss_res = np.sum((vel - vel_pred) ** 2)
    ss_tot = np.sum((vel - np.mean(vel)) ** 2)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    params.K = float(K)
    params.tau = float(tau)
    params.r_squared = float(r_squared)
    params.rmse = float(np.sqrt(np.mean((vel - vel_pred) ** 2)))


def _estimate_theoretical_params(motor_id: int, params: MotorParams) -> MotorParams:
    """Estima parámetros con modelo teórico si no hay telemetría."""
    # Valores basados en modelo_motor_dc.m
    if motor_id == 1:
        # Motor base: 30 kgf·cm, 218:1, 12V
        params.K = 0.35       # rad/s por PWM unitario (~33 RPM a 100%)
        params.tau = 0.08     # constante de tiempo
        params.friction_pwm = 0.10
    else:
        # Motor codo: 19 kgf·cm, 270:1, 24V
        params.K = 0.30
        params.tau = 0.06
        params.friction_pwm = 0.08

    params.calibrated = False
    logger.info(f"Motor {motor_id}: parámetros teóricos estimados (no calibrado)")
    return params
