# Firmware ESP32 para Robot 2R

Firmware en Rust (`no_std`) para la ESP32 que controla dos motores DC siguiendo trayectorias recibidas desde el PC.

## Arquitectura

```
┌─────────────────────────────────────────────┐
│  main loop (1 kHz)                          │
│  ├── Lee encoders (QEI)                     │
│  ├── Calcula PID                             │
│  ├── Actualiza PWM                           │
│  └── Cada 10 ms: envía telemetría            │
│                                              │
│  Tarea WiFi (background)                     │
│  ├── Acepta conexiones TCP                   │
│  ├── Recibe chunks de trayectoria            │
│  └── Responde a UDP broadcast discovery       │
└─────────────────────────────────────────────┘
```

## Compilación y flasheo

```bash
# Instalar herramientas ESP32 Rust
cargo install espup espflash
espup install

# Compilar
cd backend/esp32_firmware
cargo build --release

# Flashear
espflash flash --monitor target/xtensa-esp32-espidf/release/esp32_2r_firmware
```

## Configuración de pines

| Función       | Pin ESP32 | Notas                  |
|---------------|-----------|------------------------|
| Motor 1 PWM   | GPIO 25   | Canal LEDC 0           |
| Motor 1 DIR   | GPIO 26   | HIGH = adelante        |
| Motor 2 PWM   | GPIO 27   | Canal LEDC 1           |
| Motor 2 DIR   | GPIO 14   | HIGH = adelante        |
| Encoder 1 A   | GPIO 34   | QEI entrada A          |
| Encoder 1 B   | GPIO 35   | QEI entrada B          |
| Encoder 2 A   | GPIO 36   | QEI entrada A          |
| Encoder 2 B   | GPIO 39   | QEI entrada B          |
| UART0 (USB)   | TX0/RX0   | Comunicación con PC    |

## Protocolo

Mismo protocolo binario que `rust_comm` (ver `backend/rust_comm/src/protocol.rs`):
- Sync word: `0xAA 0x55`
- Trama con CRC32
- Tipos: NetworkConfig, TrajectoryChunk, StartExecution, Telemetry, etc.
