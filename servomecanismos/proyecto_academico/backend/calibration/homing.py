"""
Calibración de posición Home y cero del robot 2R.

El proceso de homing establece:
1. La posición de "cero" de cada articulación (offset del encoder).
2. Los límites angulares de cada motor.
3. La capacidad de recordar la posición entre ciclos de encendido.

Flujo:
    1. El usuario posiciona manualmente el robot en su posición home conocida
       (base: θ₁=175°, codo: θ₂=-128° replegado debajo del trébol).
    2. Se registra el offset del encoder para cada motor.
    3. Se buscan los límites angulares moviendo con PWM mínimo.
    4. Los offsets se guardan en NVS de la ESP32 para persistencia.

Si la ESP32 no puede recordar la posición (NVS no disponible o corrupto),
se fuerza una recalibración al inicio.
"""

from dataclasses import dataclass, field
from typing import Optional, Callable
import numpy as np
import logging

logger = logging.getLogger(__name__)


# ── Posición Home por defecto ─────────────────────────
# Robot replegado: por debajo de la línea media del trébol
# y completamente a la izquierda del mismo.
# l1 = 195 mm, l2 = 260 mm centro a centro.
HOME_TH1_DEG = 190.0    # Base: apuntando a la izquierda (x < 0.075m)
HOME_TH2_DEG = -130.0   # Codo: plegado hacia abajo (y < 0.200m)

# Tolerancia para considerar posición home alcanzada [rad]
HOME_TOLERANCE_RAD = 0.05  # ~2.86°


@dataclass
class HomeResult:
    """Resultado del proceso de homing."""
    success: bool = False
    motor1_offset: float = 0.0     # Offset del encoder motor base [rad]
    motor2_offset: float = 0.0     # Offset del encoder motor codo [rad]
    motor1_angle_min: float = -np.pi
    motor1_angle_max: float = np.pi
    motor2_angle_min: float = -np.pi
    motor2_angle_max: float = np.pi
    needs_recalibration: bool = False
    message: str = ""


class HomeCalibration:
    """
    Orquestador del proceso de calibración de posición home.

    Uso:
        home = HomeCalibration(comm_mgr)
        home.start_base_calibration()    # Calibra motor base
        home.start_elbow_calibration()   # Calibra motor codo (con base en home)
    """

    def __init__(self, comm_manager=None):
        self._comm = comm_manager
        self._state = "idle"
        self._progress_callback: Optional[Callable[[str, float], None]] = None

    @property
    def home_th1_rad(self) -> float:
        return np.deg2rad(HOME_TH1_DEG)

    @property
    def home_th2_rad(self) -> float:
        return np.deg2rad(HOME_TH2_DEG)

    def set_progress_callback(self, cb: Callable[[str, float], None]):
        """Establece callback para reportar progreso a la GUI."""
        self._progress_callback = cb

    def _report(self, status: str, progress: float = 0.0):
        if self._progress_callback:
            self._progress_callback(status, progress)

    # ── Fase 1: Calibración del motor base ────────────

    def calibrate_base_zero(
        self,
        current_theta1: float,
        known_home_theta1: Optional[float] = None,
    ) -> float:
        """
        Calcula el offset del encoder de la base.

        El usuario coloca el robot en la posición home conocida
        (o indica cuál es el ángulo actual). Se calcula:
            offset = current_encoder_value - known_home_angle

        Args:
            current_theta1: Lectura actual del encoder [rad].
            known_home_theta1: Ángulo home conocido (default: HOME_TH1_DEG).

        Returns:
            Offset del encoder [rad].
        """
        home_th1 = known_home_theta1 if known_home_theta1 is not None else self.home_th1_rad
        offset = current_theta1 - home_th1
        logger.info(f"Base: offset calculado = {np.rad2deg(offset):.2f}° "
                    f"(encoder={np.rad2deg(current_theta1):.2f}°, home={np.rad2deg(home_th1):.2f}°)")
        return offset

    # ── Fase 2: Calibración del motor codo ─────────────
    # (requiere que la base esté en posición home)

    def calibrate_elbow_zero(
        self,
        current_theta2: float,
        known_home_theta2: Optional[float] = None,
    ) -> float:
        """
        Calcula el offset del encoder del codo.

        PRECONDICIÓN: La base debe estar en su posición home.

        Args:
            current_theta2: Lectura actual del encoder [rad].
            known_home_theta2: Ángulo home conocido (default: HOME_TH2_DEG).

        Returns:
            Offset del encoder [rad].
        """
        home_th2 = known_home_theta2 if known_home_theta2 is not None else self.home_th2_rad
        offset = current_theta2 - home_th2
        logger.info(f"Codo: offset calculado = {np.rad2deg(offset):.2f}° "
                    f"(encoder={np.rad2deg(current_theta2):.2f}°, home={np.rad2deg(home_th2):.2f}°)")
        return offset

    # ── Verificación de home ──────────────────────────

    def is_at_home(self, theta1: float, theta2: float) -> bool:
        """Verifica si el robot está en la posición home."""
        d1 = abs(theta1 - self.home_th1_rad)
        d2 = abs(theta2 - self.home_th2_rad)
        # Desenvolver ángulos
        d1 = min(d1, 2 * np.pi - d1)
        d2 = min(d2, 2 * np.pi - d2)
        return d1 < HOME_TOLERANCE_RAD and d2 < HOME_TOLERANCE_RAD

    # ── Persistencia ──────────────────────────────────

    def save_home_to_esp32(self, offsets: dict[str, float]) -> bool:
        """
        Guarda los offsets de home en la memoria persistente de la ESP32.

        Args:
            offsets: {'motor1_offset': float, 'motor2_offset': float}

        Returns:
            True si se guardó exitosamente.
        """
        if self._comm is None:
            logger.warning("Sin comunicación, no se puede guardar en ESP32")
            return False

        try:
            result = self._comm.save_calibration(offsets)
            logger.info("Offsets guardados en ESP32 NVS")
            return result
        except Exception as e:
            logger.error(f"Error guardando offsets en ESP32: {e}")
            return False

    def load_home_from_esp32(self) -> Optional[dict[str, float]]:
        """
        Carga los offsets de home desde la memoria persistente de la ESP32.

        Returns:
            Dict con offsets o None si no están disponibles.
        """
        if self._comm is None:
            return None

        try:
            data = self._comm.load_calibration()
            if data and data.get('motor1_offset') is not None:
                logger.info("Offsets cargados desde ESP32 NVS")
                return data
        except Exception as e:
            logger.warning(f"No se pudieron cargar offsets de ESP32: {e}")

        return None
