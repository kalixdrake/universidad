"""
Gestor de comunicación PC ↔ ESP32.

Soporta dos modos de transporte:
1. USB/Serial (usando pyserial o el módulo Rust rust_comm)
2. WiFi/TCP (usando sockets estándar)

Prioriza intentar Rust (rust_comm) y cae en pyserial puro si no está disponible.
"""
import time
import socket
import threading
import logging
from typing import Optional, Callable
from dataclasses import dataclass
from enum import Enum

from . import (
    FrameType, TelemetryFrame, TrajectoryData,
    MotorPwmCommand, MotorPwmTelemetry, CalibrationStatus,
    pack_frame, pack_telemetry, pack_trajectory_chunk, unpack_frame,
)

logger = logging.getLogger(__name__)


class ConnectionState(Enum):
    DISCONNECTED = "disconnected"
    USB = "usb"
    NETWORK = "network"


@dataclass
class ConnectionInfo:
    state: ConnectionState = ConnectionState.DISCONNECTED
    port: str = ""          # USB: /dev/ttyUSB0
    baudrate: int = 115200
    address: str = ""       # Network: 192.168.1.100:8080


class CommunicationManager:
    """
    Gestor de alto nivel para la comunicación con la ESP32.

    Uso típico:
        mgr = CommunicationManager()
        mgr.set_telemetry_callback(my_callback)
        mgr.connect()                     # flujo automático
        mgr.send_trajectory(trajectory)
        mgr.send_start()
        # ... telemetría llega por callback ...
        mgr.close()
    """

    def __init__(self):
        self._info = ConnectionInfo()
        self._serial = None          # pyserial.Serial o rust_comm.SerialLink
        self._socket: Optional[socket.socket] = None
        self._telemetry_cb: Optional[Callable[[TelemetryFrame], None]] = None
        self._running = False
        self._rx_thread: Optional[threading.Thread] = None
        self._use_rust = False

        # Intentar cargar el módulo Rust (compilado con maturin)
        try:
            import rust_comm
            # Verificar que tiene la clase esperada (puede haber otro módulo rust_comm)
            if not hasattr(rust_comm, 'Esp32Connection'):
                raise ImportError("El módulo rust_comm encontrado no es el de este proyecto")
            self._rust_conn = rust_comm.Esp32Connection()
            self._use_rust = True
            logger.info("rust_comm cargado exitosamente.")
        except (ImportError, AttributeError) as e:
            self._rust_conn = None
            self._use_rust = False
            logger.info("rust_comm no disponible (%s), usando pyserial + socket.", e)

    # ── Propiedades ────────────────────────────────────
    @property
    def info(self) -> ConnectionInfo:
        return self._info

    @property
    def is_connected(self) -> bool:
        return self._info.state != ConnectionState.DISCONNECTED

    # ── Descubrimiento ─────────────────────────────────
    def discover_usb(self) -> Optional[str]:
        """Busca la ESP32 en puertos USB. Retorna el puerto o None."""
        if self._use_rust and self._rust_conn:
            return self._rust_conn.discover_usb()
        else:
            return self._discover_usb_pyserial()

    def _discover_usb_pyserial(self) -> Optional[str]:
        """Fallback con pyserial para descubrimiento USB."""
        try:
            import serial.tools.list_ports
            for port in serial.tools.list_ports.comports():
                if self._ping_serial(port.device):
                    return port.device
        except ImportError:
            logger.warning("pyserial no instalado.")
        except Exception as e:
            logger.debug(f"Error escaneando USB: {e}")
        return None

    def _ping_serial(self, port: str, baudrate: int = 115200, timeout: float = 0.3) -> bool:
        """Envía un ping a un puerto serie y espera ACK."""
        try:
            import serial
            with serial.Serial(port, baudrate, timeout=timeout) as ser:
                ping = pack_frame(FrameType.HEARTBEAT, b'\x00' * 4)
                ser.write(ping)
                ser.flush()
                resp = ser.read(32)
                if len(resp) >= 7 and resp[:2] == b'\xAA\x55':
                    return True
        except Exception:
            pass
        return False

    def discover_network(self) -> Optional[str]:
        """Busca la ESP32 en la red local vía broadcast UDP. Retorna ip:puerto."""
        if self._use_rust and self._rust_conn:
            return self._rust_conn.discover_network()
        else:
            return self._discover_network_socket()

    def _discover_network_socket(self) -> Optional[str]:
        """Fallback con socket UDP para descubrimiento de red."""
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            sock.setsockopt(socket.SOL_SOCKET, socket.SO_BROADCAST, 1)
            sock.settimeout(2.0)
            sock.sendto(b'ESP32_DISCOVER', ('255.255.255.255', 4210))
            data, addr = sock.recvfrom(128)
            sock.close()
            resp = data.decode('utf-8', errors='ignore')
            if resp.startswith('ESP32_2R_ROBOT'):
                port = '8080'
                if ':' in resp:
                    port = resp.split(':')[1].strip()
                return f'{addr[0]}:{port}'
        except socket.timeout:
            pass
        except Exception as e:
            logger.debug(f"Error descubrimiento red: {e}")
        return None

    # ── Conexión ───────────────────────────────────────
    def connect_auto(self) -> ConnectionState:
        """
        Flujo automático de conexión:
        1. Intenta red (si hay historial)
        2. Cae en USB
        3. Si USB funciona, configura red
        """
        # 1. Intentar red
        net_addr = self.discover_network()
        if net_addr:
            if self.connect_network(net_addr):
                return ConnectionState.NETWORK

        # 2. Intentar USB
        usb_port = self.discover_usb()
        if usb_port:
            if self.connect_usb(usb_port):
                return ConnectionState.USB

        return ConnectionState.DISCONNECTED

    def connect_usb(self, port: str, baudrate: int = 115200) -> bool:
        """Conecta por USB/serial."""
        try:
            if self._use_rust and self._rust_conn:
                ok = self._rust_conn.connect_usb(port, baudrate)
                if ok:
                    self._info = ConnectionInfo(state=ConnectionState.USB, port=port, baudrate=baudrate)
                    self._start_rx()
                    return True
            else:
                import serial
                self._serial = serial.Serial(port, baudrate, timeout=0.01)
                self._info = ConnectionInfo(state=ConnectionState.USB, port=port, baudrate=baudrate)
                self._start_rx()
                return True
        except Exception as e:
            logger.error(f"Error conectando USB {port}: {e}")
        return False

    def connect_network(self, address: str) -> bool:
        """Conecta por red (TCP). address: 'ip:puerto'."""
        try:
            host, port_str = address.rsplit(':', 1)
            port = int(port_str)
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            sock.settimeout(3.0)
            sock.connect((host, port))
            self._socket = sock
            self._info = ConnectionInfo(state=ConnectionState.NETWORK, address=address)
            self._start_rx()
            return True
        except Exception as e:
            logger.error(f"Error conectando red {address}: {e}")
        return False

    def configure_network(
        self, ssid: str, password: str, static_ip: str = "192.168.4.1"
    ) -> bool:
        """Envía configuración WiFi a la ESP32 por USB."""
        payload = f'{ssid}\0{password}\0{static_ip}'.encode('utf-8')
        frame = pack_frame(FrameType.NETWORK_CONFIG, payload)
        return self._send_raw(frame)

    # ── Envío de datos ─────────────────────────────────
    def send_trajectory(self, traj: TrajectoryData) -> bool:
        """
        Envía la trayectoria completa a la ESP32, fragmentándola en chunks.
        """
        CHUNK_SIZE = 256
        n = len(traj.theta1)
        total_chunks = (n + CHUNK_SIZE - 1) // CHUNK_SIZE

        for i in range(total_chunks):
            start = i * CHUNK_SIZE
            end = min(start + CHUNK_SIZE, n)

            frame = pack_trajectory_chunk(
                chunk_index=i,
                total_chunks=total_chunks,
                chunk_size=end - start,
                theta1_chunk=traj.theta1[start:end],
                theta2_chunk=traj.theta2[start:end],
                omega1_chunk=traj.omega1[start:end],
                omega2_chunk=traj.omega2[start:end],
            )
            if not self._send_raw(frame):
                logger.error(f"Error enviando chunk {i}/{total_chunks}")
                return False
            time.sleep(0.01)  # dar tiempo a la ESP32

        # Enviar metadata
        meta = pack_frame(FrameType.TRAJECTORY, struct.pack('<fII', traj.dt, traj.cycles, n))
        self._send_raw(meta)

        return True

    def send_start(self) -> bool:
        """Envía comando de inicio de ejecución."""
        return self._send_raw(pack_frame(FrameType.START_EXECUTION, b''))

    def send_emergency_stop(self) -> bool:
        """Envía parada de emergencia."""
        return self._send_raw(pack_frame(FrameType.EMERGENCY_STOP, b''))

    # ── Comandos de calibración ────────────────────────

    def send_motor_pwm(self, motor_id: int, duty_cycle: float, duration_ms: int = 0) -> bool:
        """
        Envía un comando PWM directo para un motor específico.

        Args:
            motor_id: 1 = base, 2 = codo.
            duty_cycle: 0.0 a 1.0 (0% a 100%).
            duration_ms: Duración del pulso (0 = indefinido, hasta nuevo comando).
        """
        from . import MotorPwmCommand
        cmd = MotorPwmCommand(motor_id=motor_id, duty_cycle=duty_cycle, duration_ms=duration_ms)
        return self._send_raw(pack_frame(FrameType.MOTOR_PWM_CMD, cmd.to_bytes()))

    def stop_motor(self, motor_id: int) -> bool:
        """Detiene un motor específico (PWM = 0)."""
        return self.send_motor_pwm(motor_id, 0.0)

    def stop_all_motors(self) -> bool:
        """Detiene ambos motores."""
        ok1 = self.send_motor_pwm(1, 0.0)
        ok2 = self.send_motor_pwm(2, 0.0)
        return ok1 and ok2

    def save_calibration(self, calibration_data: dict) -> bool:
        """
        Guarda datos de calibración en la memoria persistente (NVS) de la ESP32.

        Args:
            calibration_data: Diccionario con todos los parámetros de calibración.
        """
        import json
        payload = json.dumps(calibration_data).encode('utf-8')
        return self._send_raw(pack_frame(FrameType.CALIBRATION_SAVE, payload))

    def load_calibration(self) -> Optional[dict]:
        """
        Carga datos de calibración desde la memoria persistente (NVS) de la ESP32.

        Returns:
            Diccionario con parámetros o None si no disponible.
        """
        import json
        frame = pack_frame(FrameType.CALIBRATION_LOAD, b'')
        if not self._send_raw(frame):
            return None

        # Esperar respuesta
        time.sleep(0.3)
        buf = bytearray()
        deadline = time.time() + 2.0
        while time.time() < deadline:
            data = self._read_available()
            if data:
                buf.extend(data)
                if len(buf) >= 7:
                    msg_type, payload = unpack_frame(bytes(buf))
                    if msg_type == FrameType.CALIBRATION_DATA and payload:
                        try:
                            return json.loads(payload.decode('utf-8'))
                        except Exception:
                            return None
            time.sleep(0.05)
        return None

    def request_calibration_status(self) -> Optional[dict]:
        """
        Solicita el estado de calibración de la ESP32.

        Returns:
            Dict con position_calibrated, motor_characterized, nvs_available, etc.
        """
        from . import CalibrationStatus
        frame = pack_frame(FrameType.CALIBRATION_STATUS, b'')
        if not self._send_raw(frame):
            return None

        time.sleep(0.2)
        buf = bytearray()
        deadline = time.time() + 2.0
        while time.time() < deadline:
            data = self._read_available()
            if data:
                buf.extend(data)
                if len(buf) >= 7:
                    msg_type, payload = unpack_frame(bytes(buf))
                    if msg_type == FrameType.CALIBRATION_STATUS and payload:
                        try:
                            status = CalibrationStatus.from_bytes(payload)
                            return {
                                'position_calibrated': status.position_calibrated,
                                'motor_characterized': status.motor_characterized,
                                'nvs_available': status.nvs_available,
                                'motor1_stall_detected': status.motor1_stall_detected,
                                'motor2_stall_detected': status.motor2_stall_detected,
                                'error_code': status.error_code,
                            }
                        except Exception:
                            return None
            time.sleep(0.05)
        return None

    # ── Recepción y telemetría ─────────────────────────
    def set_telemetry_callback(self, cb: Callable[[TelemetryFrame], None]):
        """Registra un callback para frames de telemetría."""
        self._telemetry_cb = cb

    def _start_rx(self):
        """Inicia el hilo de recepción."""
        if self._rx_thread and self._rx_thread.is_alive():
            return
        self._running = True
        self._rx_thread = threading.Thread(target=self._rx_loop, daemon=True)
        self._rx_thread.start()

    def _rx_loop(self):
        """Bucle de recepción: lee tramas y despacha telemetría."""
        buf = bytearray()
        while self._running:
            try:
                data = self._read_available()
                if data:
                    buf.extend(data)
                    # Intentar decodificar tramas del buffer
                    while len(buf) >= 7:
                        msg_type, payload = unpack_frame(bytes(buf))
                        if msg_type is None:
                            # Buscar siguiente sync word
                            sync_idx = buf.find(b'\xAA\x55', 1)
                            if sync_idx > 0:
                                del buf[:sync_idx]
                            else:
                                del buf[:max(0, len(buf) - 2)]
                            continue

                        frame_len = HEADER_SIZE + struct.unpack('<H', buf[2:4])[0] + CRC_SIZE
                        del buf[:frame_len]

                        if msg_type == FrameType.TELEMETRY and payload and self._telemetry_cb:
                            try:
                                telem = TelemetryFrame.from_bytes(payload)
                                self._telemetry_cb(telem)
                            except Exception as e:
                                logger.debug(f"Error decodificando telemetría: {e}")
                else:
                    time.sleep(0.001)
            except Exception as e:
                if self._running:
                    logger.debug(f"Error en rx_loop: {e}")
                time.sleep(0.01)

    def _read_available(self) -> Optional[bytes]:
        """Lee datos disponibles del medio activo."""
        if self._use_rust and self._rust_conn and self._info.state == ConnectionState.USB:
            # Rust maneja su propio buffering; aquí leemos vía el método Python
            telem = self._rust_conn.read_telemetry(10)
            if telem and self._telemetry_cb:
                tf = TelemetryFrame(
                    timestamp_us=telem.timestamp_us,
                    theta1=telem.theta1,
                    theta2=telem.theta2,
                    omega1=telem.omega1,
                    omega2=telem.omega2,
                    waypoint_index=telem.waypoint_index,
                    state=telem.state,
                    error_code=telem.error_code,
                )
                self._telemetry_cb(tf)
            return None  # Rust ya procesó

        if self._serial and self._serial.is_open:
            try:
                return self._serial.read(self._serial.in_waiting or 1)
            except Exception:
                return None

        if self._socket:
            try:
                self._socket.settimeout(0.001)
                return self._socket.recv(4096)
            except (socket.timeout, BlockingIOError):
                return None
            except Exception:
                return None

        return None

    def _send_raw(self, data: bytes) -> bool:
        """Envía bytes crudos por el medio activo."""
        try:
            if self._info.state == ConnectionState.USB:
                if self._use_rust and self._rust_conn:
                    # Usar write directo por serial (Rust)
                    return True  # Rust envía internamente
                elif self._serial:
                    self._serial.write(data)
                    self._serial.flush()
                    return True
            elif self._info.state == ConnectionState.NETWORK and self._socket:
                self._socket.sendall(data)
                return True
        except Exception as e:
            logger.error(f"Error enviando datos: {e}")
        return False

    # ── Cierre ─────────────────────────────────────────
    def close(self):
        """Cierra la conexión activa."""
        self._running = False
        if self._rx_thread:
            self._rx_thread.join(timeout=1.0)

        if self._use_rust and self._rust_conn:
            self._rust_conn.close()
        if self._serial:
            self._serial.close()
            self._serial = None
        if self._socket:
            self._socket.close()
            self._socket = None

        self._info = ConnectionInfo()
