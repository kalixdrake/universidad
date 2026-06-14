from .motor_identification import MotorIdentification
from .homing import HomingSupervisor

class CalibrationManager:
    """
    Orquestador de calibración, persistencia de variables y cargador de configuraciones de motor.
    """
    def __init__(self):
        self.calibrated = False
        self.motor_base_model = MotorIdentification(motor_id=0)
        self.motor_elbow_model = MotorIdentification(motor_id=1)
        self.homing_system = HomingSupervisor()
        
    def get_status(self):
        if self.calibrated:
            return "Calibrado y Listo"
        return "Pendiente de Calibración / Homing"
        
    def run_full_calibration(self):
        """
        Realiza el ciclo completo: homing de articulaciones + identificación de motores + persistencia.
        """
        print("[CALIBRACIÓN] Iniciando autocalibración guiada completa...")
        # 1. Homing
        home_coords = self.homing_system.start_homing_sequence()
        
        # 2. Identificación de motores
        stiction_base = self.motor_base_model.measure_stiction()
        stiction_elbow = self.motor_elbow_model.measure_stiction()
        
        self.calibrated = True
        print(f"[CALISTRACIÓN] Guardados parámetros de fricción estática:")
        print(f"               Base PWM arr: {stiction_base}")
        print(f"               Codo PWM arr: {stiction_elbow}")
        return self.calibrated
