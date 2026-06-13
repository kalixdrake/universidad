"""
Protocolo de comunicación PC ↔ ESP32 para robot 2R.

Define las estructuras de datos que viajan entre el backend y la ESP32,
independientemente del medio de transporte (USB/serial o WiFi/TCP).
"""
from dataclasses import dataclass, field
from enum import IntEnum
from typing import Optional
import struct
import zlib


class FrameType(IntEnum):
    NETWORK_CONFIG = 0x01
    TRAJECTORY = 0x02
    TRAJECTORY_CHUNK = 0x03
    START_EXECUTION = 0x10
    EMERGENCY_STOP = 0x11
    ACK = 0x20
    TELEMETRY = 0x30
    HEARTBEAT = 0x40
    # ── Comandos de calibración ───────────────────────
    CALIBRATION_CMD = 0x50     # Comando de calibración genérico
    CALIBRATION_DATA = 0x51    # Respuesta con datos de calibración
    CALIBRATION_SAVE = 0x52    # Guardar calibración en NVS
    CALIBRATION_LOAD = 0x53    # Cargar calibración de NVS
    CALIBRATION_STATUS = 0x54  # Estado de calibración
    MOTOR_PWM_CMD = 0x60       # Control PWM directo por motor
    MOTOR_PWM_TELEM = 0x61     # Telemetría extendida con PWM actual


SYNC_WORD = b'\xAA\x55'
HEADER_SIZE = 5   # sync(2) + length(2) + type(1)
CRC_SIZE = 4
MAX_PAYLOAD = 16 * 1024


@dataclass
class TelemetryFrame:
    """Telemetría en tiempo real desde la ESP32."""
    timestamp_us: int = 0
    theta1: float = 0.0
    theta2: float = 0.0
    omega1: float = 0.0
    omega2: float = 0.0
    waypoint_index: int = 0
    state: int = 0       # 0=idle, 1=running, 2=error, 3=completed
    error_code: int = 0

    def to_bytes(self) -> bytes:
        return struct.pack(
            '<QddddIBB',
            self.timestamp_us,
            self.theta1, self.theta2,
            self.omega1, self.omega2,
            self.waypoint_index,
            self.state, self.error_code,
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> 'TelemetryFrame':
        vals = struct.unpack('<QddddIBB', data[:50])
        return cls(
            timestamp_us=vals[0],
            theta1=vals[1], theta2=vals[2],
            omega1=vals[3], omega2=vals[4],
            waypoint_index=vals[5],
            state=vals[6], error_code=vals[7],
        )


@dataclass
class TrajectoryData:
    """Datos completos de trayectoria (se fragmentan para envío)."""
    theta1: list[float] = field(default_factory=list)
    theta2: list[float] = field(default_factory=list)
    omega1: list[float] = field(default_factory=list)
    omega2: list[float] = field(default_factory=list)
    dt: float = 0.01
    cycles: int = 1


# ── Datos de calibración ──────────────────────────────

@dataclass
class MotorPwmCommand:
    """Comando PWM directo para un motor específico."""
    motor_id: int = 1        # 1=base, 2=codo
    duty_cycle: float = 0.0  # 0.0 a 1.0
    duration_ms: int = 0     # 0 = indefinido

    def to_bytes(self) -> bytes:
        return struct.pack('<BfI', self.motor_id, self.duty_cycle, self.duration_ms)

    @classmethod
    def from_bytes(cls, data: bytes) -> 'MotorPwmCommand':
        motor_id, duty, dur = struct.unpack('<BfI', data[:9])
        return cls(motor_id=motor_id, duty_cycle=duty, duration_ms=dur)


@dataclass
class MotorPwmTelemetry:
    """Telemetría extendida incluyendo PWM actual y encoder raw."""
    timestamp_us: int = 0
    theta1: float = 0.0
    theta2: float = 0.0
    omega1: float = 0.0
    omega2: float = 0.0
    pwm1: float = 0.0          # PWM actual motor 1
    pwm2: float = 0.0          # PWM actual motor 2
    encoder1_raw: int = 0      # Encoder raw motor 1
    encoder2_raw: int = 0      # Encoder raw motor 2

    def to_bytes(self) -> bytes:
        return struct.pack(
            '<Qddddffii',
            self.timestamp_us,
            self.theta1, self.theta2,
            self.omega1, self.omega2,
            self.pwm1, self.pwm2,
            self.encoder1_raw, self.encoder2_raw,
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> 'MotorPwmTelemetry':
        vals = struct.unpack('<Qddddffii', data[:52])
        return cls(
            timestamp_us=vals[0],
            theta1=vals[1], theta2=vals[2],
            omega1=vals[3], omega2=vals[4],
            pwm1=vals[5], pwm2=vals[6],
            encoder1_raw=vals[7], encoder2_raw=vals[8],
        )


@dataclass
class CalibrationStatus:
    """Estado de calibración reportado por la ESP32."""
    position_calibrated: bool = False
    motor_characterized: bool = False
    nvs_available: bool = False
    motor1_stall_detected: bool = False
    motor2_stall_detected: bool = False
    error_code: int = 0

    def to_bytes(self) -> bytes:
        flags = (
            (1 if self.position_calibrated else 0) |
            (2 if self.motor_characterized else 0) |
            (4 if self.nvs_available else 0) |
            (8 if self.motor1_stall_detected else 0) |
            (16 if self.motor2_stall_detected else 0)
        )
        return struct.pack('<BB', flags, self.error_code)

    @classmethod
    def from_bytes(cls, data: bytes) -> 'CalibrationStatus':
        flags, err = struct.unpack('<BB', data[:2])
        return cls(
            position_calibrated=bool(flags & 1),
            motor_characterized=bool(flags & 2),
            nvs_available=bool(flags & 4),
            motor1_stall_detected=bool(flags & 8),
            motor2_stall_detected=bool(flags & 16),
            error_code=err,
        )


# ── Empaquetado ──────────────────────────────────────

def pack_frame(msg_type: FrameType, payload: bytes) -> bytes:
    """Empaqueta una trama con sync word, longitud, tipo y CRC32."""
    length = len(payload)
    header = SYNC_WORD + struct.pack('<HB', length, msg_type.value)
    crc = zlib.crc32(header + payload) & 0xFFFFFFFF
    return header + payload + struct.pack('<I', crc)


def unpack_frame(data: bytes) -> tuple[Optional[FrameType], Optional[bytes]]:
    """Desempaqueta una trama verificando CRC. Retorna (tipo, payload) o (None, None)."""
    if len(data) < HEADER_SIZE + CRC_SIZE:
        return None, None
    if data[:2] != SYNC_WORD:
        return None, None

    length = struct.unpack('<H', data[2:4])[0]
    msg_type_byte = data[4]
    payload_end = HEADER_SIZE + length

    if len(data) < payload_end + CRC_SIZE:
        return None, None

    payload = data[HEADER_SIZE:payload_end]
    expected_crc = struct.unpack('<I', data[payload_end:payload_end + CRC_SIZE])[0]
    actual_crc = zlib.crc32(data[:payload_end]) & 0xFFFFFFFF

    if actual_crc != expected_crc:
        return None, None

    try:
        msg_type = FrameType(msg_type_byte)
    except ValueError:
        return None, None

    return msg_type, payload


def pack_telemetry(telem: TelemetryFrame) -> bytes:
    """Empaqueta una trama de telemetría."""
    return pack_frame(FrameType.TELEMETRY, telem.to_bytes())


def pack_trajectory_chunk(
    chunk_index: int,
    total_chunks: int,
    chunk_size: int,
    theta1_chunk: list[float],
    theta2_chunk: list[float],
    omega1_chunk: list[float],
    omega2_chunk: list[float],
) -> bytes:
    """Empaqueta un chunk de trayectoria."""
    n = len(theta1_chunk)
    fmt = f'<II{3*n}d'
    # Pack theta1, theta2, omega1, omega2 interleaved
    payload = bytearray()
    payload.extend(struct.pack('<III', chunk_index, total_chunks, n))
    for i in range(n):
        payload.extend(struct.pack('<dddd',
            theta1_chunk[i], theta2_chunk[i],
            omega1_chunk[i], omega2_chunk[i]))
    return pack_frame(FrameType.TRAJECTORY_CHUNK, bytes(payload))
