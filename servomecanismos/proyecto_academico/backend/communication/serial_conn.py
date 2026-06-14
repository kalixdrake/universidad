import serial
import time

class SerialConnection:
    """
    Gestiona la conexión serie (UART) entre la PC y la ESP32.
    """
    def __init__(self, port: str = "/dev/ttyUSB0", baudrate: int = 115200, timeout: float = 1.0):
        self.port = port
        self.baudrate = baudrate
        self.timeout = timeout
        self.ser = None
        
    def connect(self):
        try:
            self.ser = serial.Serial(self.port, self.baudrate, timeout=self.timeout)
            # Permite el reinicio de la ESP32 al conectar
            time.sleep(1.0)
            return True
        except Exception as e:
            print(f"[ERROR] No se pudo abrir el puerto serie {self.port}: {e}")
            return False
            
    def send_command(self, cmd_type: str, data: bytes):
        """
        Envía un comando serial empaquetado a la ESP32.
        """
        if self.ser is None or not self.ser.is_open:
            return False
            
        # Esqueleto de trama UART (Formato: [Cabecera][Tipo][Largo][Datos][CRC32][Cola])
        # Esto se complementará luego de definir el protocolo final
        header = b"\xAA\x55"
        payload = header + cmd_type.encode() + len(data).to_bytes(2, "big") + data
        
        self.ser.write(payload)
        return True
        
    def read_telemetry(self):
        """
        Lee datos de telemetría de retorno enviados por la ESP32.
        """
        if self.ser is None or not self.ser.is_open:
            return None
            
        if self.ser.in_watermark > 0 or self.ser.in_waiting > 0:
            return self.ser.readline()
        return None
        
    def disconnect(self):
        if self.ser and self.ser.is_open:
            self.ser.close()
            print("[INFO] Conexión serie cerrada.")
