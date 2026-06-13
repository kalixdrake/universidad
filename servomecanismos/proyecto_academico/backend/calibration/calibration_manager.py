"""
Gestor principal de calibración del robot 2R.

Orquesta el flujo completo de calibración:
    1. Verificar si ya existe calibración guardada (ESP32 NVS).
    2. Si no existe o está corrupta → forzar recalibración completa.
    3. Si existe → ofrecer opción de recalibrar posición solamente.
    
La calibración de posición (home/cero/límites) se puede repetir.
La calibración de curvas motor (FT, caracterización) se hace una sola vez.
"""

from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Callable
import json
import os
import logging

from .motor_identification import (
    MotorParams,
    identify_friction_torque,
    find_angle_limits,
    characterize_motor_curve,
)
from .homing import HomeCalibration, HomeResult, HOME_TH1_DEG, HOME_TH2_DEG

logger = logging.getLogger(__name__)


class CalibrationState(Enum):
    """Estado global de la calibración."""
    NOT_CALIBRATED = auto()       # Nunca calibrado
    NEEDS_RECALIBRATION = auto()  # Datos corruptos/ausentes
    POSITION_ONLY = auto()        # Solo posición calibrada
    FULLY_CALIBRATED = auto()     # Todo calibrado
    IN_PROGRESS = auto()          # Calibración en curso
    ERROR = auto()                # Error durante calibración


class CalibrationStep(Enum):
    """Pasos del asistente de calibración."""
    WELCOME = auto()              # Pantalla de bienvenida
    CHECK_PREVIOUS = auto()       # Verificar calibración previa
    BASE_FRICTION = auto()        # PWM rozamiento base
    BASE_LIMITS = auto()          # Límites angulares base
    BASE_ZERO = auto()             # Cero de la base
    ELBOW_FRICTION = auto()       # PWM rozamiento codo
    ELBOW_LIMITS = auto()         # Límites angulares codo
    ELBOW_ZERO = auto()           # Cero del codo
    MOTOR_CHARACTERIZATION = auto()  # Caracterización de curvas
    SAVE = auto()                 # Guardar en ESP32
    COMPLETE = auto()             # Finalizado


@dataclass
class CalibrationData:
    """Datos completos de calibración del robot."""
    motor1: MotorParams = field(default_factory=MotorParams)
    motor2: MotorParams = field(default_factory=MotorParams)

    # Home offsets (para restaurar posición cero)
    home_th1_offset: float = 0.0
    home_th2_offset: float = 0.0

    # Timestamps
    position_calibrated_at: str = ""
    motor_characterized_at: str = ""

    # Versión del formato de calibración
    version: int = 1

    def to_dict(self) -> dict:
        return {
            'version': self.version,
            'motor1': self.motor1.to_dict(),
            'motor2': self.motor2.to_dict(),
            'home_th1_offset': self.home_th1_offset,
            'home_th2_offset': self.home_th2_offset,
            'position_calibrated_at': self.position_calibrated_at,
            'motor_characterized_at': self.motor_characterized_at,
        }

    @classmethod
    def from_dict(cls, data: dict) -> 'CalibrationData':
        return cls(
            version=data.get('version', 1),
            motor1=MotorParams.from_dict(data.get('motor1', {})),
            motor2=MotorParams.from_dict(data.get('motor2', {})),
            home_th1_offset=data.get('home_th1_offset', 0.0),
            home_th2_offset=data.get('home_th2_offset', 0.0),
            position_calibrated_at=data.get('position_calibrated_at', ''),
            motor_characterized_at=data.get('motor_characterized_at', ''),
        )

    @property
    def is_position_calibrated(self) -> bool:
        return bool(self.position_calibrated_at)

    @property
    def is_motor_characterized(self) -> bool:
        return bool(self.motor_characterized_at) and self.motor1.calibrated and self.motor2.calibrated

    @property
    def is_fully_calibrated(self) -> bool:
        return self.is_position_calibrated and self.is_motor_characterized


class CalibrationManager:
    """
    Gestor principal de calibración.

    Uso:
        mgr = CalibrationManager(comm_manager)
        state = mgr.check_state()  # NOT_CALIBRATED, FULLY_CALIBRATED, etc.
        mgr.run_full_calibration(progress_callback, telemetry_callback)
    """

    # Archivo local de respaldo
    LOCAL_CALIBRATION_FILE = "calibration_data.json"

    def __init__(self, comm_manager=None):
        self._comm = comm_manager
        self._data = CalibrationData()
        self._state = CalibrationState.NOT_CALIBRATED
        self._current_step = CalibrationStep.WELCOME
        self._home = HomeCalibration(comm_manager)

        # Callbacks
        self._step_callback: Optional[Callable[[CalibrationStep, str, float], None]] = None
        self._state_callback: Optional[Callable[[CalibrationState], None]] = None

    # ── Callbacks ─────────────────────────────────────

    def set_step_callback(self, cb: Callable[[CalibrationStep, str, float], None]):
        """Callback para reportar progreso de pasos a la GUI."""
        self._step_callback = cb

    def set_state_callback(self, cb: Callable[[CalibrationState], None]):
        """Callback para reportar cambios de estado."""
        self._state_callback = cb

    def _report_step(self, step: CalibrationStep, message: str, progress: float = 0.0):
        if self._step_callback:
            self._step_callback(step, message, progress)

    def _set_state(self, state: CalibrationState):
        self._state = state
        if self._state_callback:
            self._state_callback(state)

    # ── Verificación de estado ────────────────────────

    def check_state(self) -> CalibrationState:
        """
        Verifica el estado de calibración:
        1. Intenta cargar de ESP32 NVS.
        2. Fallback a archivo local.
        3. Si nada existe → NOT_CALIBRATED.
        """
        # Intentar cargar de ESP32
        if self._comm:
            try:
                esp32_data = self._comm.load_calibration()
                if esp32_data:
                    self._data = CalibrationData.from_dict(esp32_data)
                    if self._data.is_fully_calibrated:
                        self._set_state(CalibrationState.FULLY_CALIBRATED)
                        return self._state
                    elif self._data.is_position_calibrated:
                        self._set_state(CalibrationState.POSITION_ONLY)
                        return self._state
            except Exception as e:
                logger.warning(f"No se pudo cargar de ESP32: {e}")

        # Fallback a archivo local
        if os.path.exists(self.LOCAL_CALIBRATION_FILE):
            try:
                with open(self.LOCAL_CALIBRATION_FILE, 'r') as f:
                    local_data = json.load(f)
                self._data = CalibrationData.from_dict(local_data)
                if self._data.is_fully_calibrated:
                    self._set_state(CalibrationState.FULLY_CALIBRATED)
                elif self._data.is_position_calibrated:
                    self._set_state(CalibrationState.POSITION_ONLY)
                else:
                    self._set_state(CalibrationState.NEEDS_RECALIBRATION)
                return self._state
            except Exception as e:
                logger.warning(f"Archivo local corrupto: {e}")

        self._set_state(CalibrationState.NOT_CALIBRATED)
        return self._state

    # ── Calibración completa ──────────────────────────

    def run_full_calibration(
        self,
        telemetry_callback: Optional[Callable[[], dict]] = None,
        user_home_th1_deg: Optional[float] = None,
        user_home_th2_deg: Optional[float] = None,
    ) -> CalibrationData:
        """
        Ejecuta el flujo completo de calibración.

        Args:
            telemetry_callback: fn() -> dict con theta1, theta2, omega1, omega2.
            user_home_th1_deg: Ángulo home base indicado por el usuario.
            user_home_th2_deg: Ángulo home codo indicado por el usuario.

        Returns:
            CalibrationData con todos los parámetros.
        """
        import datetime

        self._set_state(CalibrationState.IN_PROGRESS)
        now = datetime.datetime.now().isoformat()

        home_th1 = np.deg2rad(user_home_th1_deg) if user_home_th1_deg else self._home.home_th1_rad
        home_th2 = np.deg2rad(user_home_th2_deg) if user_home_th2_deg else self._home.home_th2_rad

        # ── Paso 1: PWM rozamiento base ───────────────
        self._current_step = CalibrationStep.BASE_FRICTION
        self._report_step(self._current_step, "Buscando torque de rozamiento del motor base...", 0.05)
        friction1 = identify_friction_torque(
            1,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, s, 0.05 + 0.1 * pwm
            ),
        )

        # ── Paso 2: Límites base ──────────────────────
        self._current_step = CalibrationStep.BASE_LIMITS
        self._report_step(self._current_step, "Buscando límites angulares del motor base...", 0.15)
        angle1_min, angle1_max = find_angle_limits(
            1, friction1,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, d, a: self._report_step(
                self._current_step, f"Base: {d}={a:.1f}°", 0.15 + 0.05 * (1 if d == 'max' else 0)
            ),
        )

        # ── Paso 3: Cero base ─────────────────────────
        self._current_step = CalibrationStep.BASE_ZERO
        self._report_step(self._current_step, "Estableciendo cero del motor base...", 0.25)
        try:
            telem = telemetry_callback() if telemetry_callback else {}
            current_th1 = float(telem.get('theta1', home_th1))
            offset1 = self._home.calibrate_base_zero(current_th1, home_th1)
        except Exception:
            offset1 = 0.0

        # ── Paso 4: PWM rozamiento codo ───────────────
        self._current_step = CalibrationStep.ELBOW_FRICTION
        self._report_step(self._current_step, "Buscando torque de rozamiento del motor codo...", 0.30)
        friction2 = identify_friction_torque(
            2,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, s, 0.30 + 0.1 * pwm
            ),
        )

        # ── Paso 5: Límites codo ──────────────────────
        self._current_step = CalibrationStep.ELBOW_LIMITS
        self._report_step(self._current_step, "Buscando límites angulares del motor codo...", 0.40)
        angle2_min, angle2_max = find_angle_limits(
            2, friction2,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, d, a: self._report_step(
                self._current_step, f"Codo: {d}={a:.1f}°", 0.40 + 0.05 * (1 if d == 'max' else 0)
            ),
        )

        # ── Paso 6: Cero codo ─────────────────────────
        self._current_step = CalibrationStep.ELBOW_ZERO
        self._report_step(self._current_step, "Estableciendo cero del motor codo...", 0.50)
        try:
            telem = telemetry_callback() if telemetry_callback else {}
            current_th2 = float(telem.get('theta2', home_th2))
            offset2 = self._home.calibrate_elbow_zero(current_th2, home_th2)
        except Exception:
            offset2 = 0.0

        # ── Paso 7: Caracterización de curvas motor ───
        self._current_step = CalibrationStep.MOTOR_CHARACTERIZATION
        self._report_step(self._current_step, "Caracterizando curvas de velocidad/PWM...", 0.55)

        params1 = characterize_motor_curve(
            1, friction1, angle1_min, angle1_max,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, f"M1: {s}", 0.55 + 0.15 * pwm
            ),
        )
        params1.zero_offset = offset1

        params2 = characterize_motor_curve(
            2, friction2, angle2_min, angle2_max,
            telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, f"M2: {s}", 0.70 + 0.15 * pwm
            ),
        )
        params2.zero_offset = offset2

        # ── Paso 8: Guardar ───────────────────────────
        self._current_step = CalibrationStep.SAVE
        self._report_step(self._current_step, "Guardando parámetros de calibración...", 0.95)

        self._data.motor1 = params1
        self._data.motor2 = params2
        self._data.home_th1_offset = offset1
        self._data.home_th2_offset = offset2
        self._data.position_calibrated_at = now
        self._data.motor_characterized_at = now

        self._save_calibration()

        self._current_step = CalibrationStep.COMPLETE
        self._report_step(self._current_step, "¡Calibración completada exitosamente!", 1.0)
        self._set_state(CalibrationState.FULLY_CALIBRATED)

        return self._data

    def run_position_only_calibration(
        self,
        telemetry_callback: Optional[Callable[[], dict]] = None,
        user_home_th1_deg: Optional[float] = None,
        user_home_th2_deg: Optional[float] = None,
    ) -> CalibrationData:
        """
        Ejecuta solo la calibración de posición (sin caracterización de motores).
        Útil para recalibraciones rápidas.
        """
        import datetime

        self._set_state(CalibrationState.IN_PROGRESS)
        now = datetime.datetime.now().isoformat()

        home_th1 = np.deg2rad(user_home_th1_deg) if user_home_th1_deg else self._home.home_th1_rad
        home_th2 = np.deg2rad(user_home_th2_deg) if user_home_th2_deg else self._home.home_th2_rad

        # PWM rozamiento base
        self._current_step = CalibrationStep.BASE_FRICTION
        self._report_step(self._current_step, "Buscando torque de rozamiento base...", 0.05)
        friction1 = identify_friction_torque(
            1, telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, s, 0.05 + 0.1 * pwm
            ),
        )

        # Límites base
        self._current_step = CalibrationStep.BASE_LIMITS
        self._report_step(self._current_step, "Buscando límites base...", 0.15)
        angle1_min, angle1_max = find_angle_limits(
            1, friction1, telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, d, a: self._report_step(
                self._current_step, f"Base {d}={a:.1f}°", 0.15 + 0.1 * (1 if d == 'max' else 0)
            ),
        )

        # Cero base
        self._current_step = CalibrationStep.BASE_ZERO
        self._report_step(self._current_step, "Estableciendo cero base...", 0.30)
        try:
            telem = telemetry_callback() if telemetry_callback else {}
            offset1 = self._home.calibrate_base_zero(
                float(telem.get('theta1', home_th1)), home_th1
            )
        except Exception:
            offset1 = 0.0

        # PWM codo
        self._current_step = CalibrationStep.ELBOW_FRICTION
        self._report_step(self._current_step, "Buscando torque de rozamiento codo...", 0.35)
        friction2 = identify_friction_torque(
            2, telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, pwm, s: self._report_step(
                self._current_step, s, 0.35 + 0.1 * pwm
            ),
        )

        # Límites codo
        self._current_step = CalibrationStep.ELBOW_LIMITS
        self._report_step(self._current_step, "Buscando límites codo...", 0.50)
        angle2_min, angle2_max = find_angle_limits(
            2, friction2, telemetry_callback=telemetry_callback,
            progress_callback=lambda mid, d, a: self._report_step(
                self._current_step, f"Codo {d}={a:.1f}°", 0.50 + 0.1 * (1 if d == 'max' else 0)
            ),
        )

        # Cero codo
        self._current_step = CalibrationStep.ELBOW_ZERO
        self._report_step(self._current_step, "Estableciendo cero codo...", 0.70)
        try:
            telem = telemetry_callback() if telemetry_callback else {}
            offset2 = self._home.calibrate_elbow_zero(
                float(telem.get('theta2', home_th2)), home_th2
            )
        except Exception:
            offset2 = 0.0

        # Guardar
        self._current_step = CalibrationStep.SAVE
        self._report_step(self._current_step, "Guardando offsets de posición...", 0.90)

        self._data.motor1.friction_pwm = friction1
        self._data.motor1.angle_min = angle1_min
        self._data.motor1.angle_max = angle1_max
        self._data.motor1.zero_offset = offset1
        self._data.motor2.friction_pwm = friction2
        self._data.motor2.angle_min = angle2_min
        self._data.motor2.angle_max = angle2_max
        self._data.motor2.zero_offset = offset2
        self._data.home_th1_offset = offset1
        self._data.home_th2_offset = offset2
        self._data.position_calibrated_at = now

        self._save_calibration()

        self._current_step = CalibrationStep.COMPLETE
        self._report_step(self._current_step, "¡Calibración de posición completada!", 1.0)

        if self._data.is_motor_characterized:
            self._set_state(CalibrationState.FULLY_CALIBRATED)
        else:
            self._set_state(CalibrationState.POSITION_ONLY)

        return self._data

    # ── Persistencia ──────────────────────────────────

    def _save_calibration(self):
        """Guarda en ESP32 NVS y en archivo local de respaldo."""
        data_dict = self._data.to_dict()

        # Guardar en ESP32
        if self._comm:
            try:
                self._comm.save_calibration(data_dict)
                logger.info("Calibración guardada en ESP32 NVS")
            except Exception as e:
                logger.warning(f"No se pudo guardar en ESP32: {e}")

        # Respaldo local
        try:
            with open(self.LOCAL_CALIBRATION_FILE, 'w') as f:
                json.dump(data_dict, f, indent=2)
            logger.info(f"Calibración guardada en {self.LOCAL_CALIBRATION_FILE}")
        except Exception as e:
            logger.error(f"Error guardando archivo local: {e}")

    def load_calibration(self) -> Optional[CalibrationData]:
        """Carga datos de calibración guardados."""
        self.check_state()
        if self._state in (CalibrationState.FULLY_CALIBRATED, CalibrationState.POSITION_ONLY):
            return self._data
        return None

    @property
    def data(self) -> CalibrationData:
        return self._data

    @property
    def state(self) -> CalibrationState:
        return self._state

    @property
    def current_step(self) -> CalibrationStep:
        return self._current_step

    @property
    def home_th1_deg(self) -> float:
        return 190.0

    @property
    def home_th2_deg(self) -> float:
        return -130.0


# Import necesario al final para evitar circular
import numpy as np
