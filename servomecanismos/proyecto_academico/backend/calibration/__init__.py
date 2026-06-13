"""
Módulo de calibración del robot 2R.

Provee:
- Identificación de torque de rozamiento (PWM mínimo para movimiento)
- Calibración de posición (cero, límites angulares)
- Caracterización de curvas motor (jerk, aceleración, velocidad vs PWM)
- Almacenamiento persistente de parámetros de control

Flujo típico:
    1. Homing: encontrar cero y límites de cada articulación
    2. Identificación: PWM mínimo (torque de rozamiento)
    3. Caracterización: curvas velocidad/PWM → función de transferencia
    4. Guardado: parámetros en almacenamiento persistente (ESP32 NVS)
"""

from .calibration_manager import (
    CalibrationManager,
    CalibrationState,
    CalibrationStep,
    MotorParams,
    CalibrationData,
)

from .motor_identification import (
    identify_friction_torque,
    find_angle_limits,
    characterize_motor_curve,
    fit_transfer_function,
)

from .homing import (
    HomeCalibration,
    HomeResult,
)

__all__ = [
    'CalibrationManager',
    'CalibrationState',
    'CalibrationStep',
    'MotorParams',
    'CalibrationData',
    'identify_friction_torque',
    'find_angle_limits',
    'characterize_motor_curve',
    'fit_transfer_function',
    'HomeCalibration',
    'HomeResult',
]
