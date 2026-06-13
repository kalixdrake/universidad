# Robot 2R — Backend de Generación de Trayectorias

Backend en Python + Rust para la generación, cálculo cinemático/dinámico y comunicación en tiempo real con una ESP32 que controla un
manipulador robótico planar de 2 grados de libertad (2R).


## Arquitectura

```
┌──────────────────────────────────────────────────────────────────┐
│                        PC (Python + Rust)                         │
│  ┌────────────┐   ┌──────────────────┐   ┌────────────────────┐  │
│  │  GUI       │   │  Trayectorias    │   │  Comunicación      │  │
│  │ customtkinter│  │  numpy / scipy   │   │  pyserial / rust   │  │
│  │  • Conexión │   │  • Trébol N-pétalos│  │  • USB (serial)   │  │
│  │  • Calibración│ │  • Bézier approach│  │  • WiFi (TCP/UDP)  │  │
│  │  • Sliders  │   │  • Cinemática inv │   │  • Protocolo CRC32 │  │
│  │  • Preview  │   │  • Dinámica inv   │   │  • Calibración    │  │
│  │  • Monitoreo │  └──────────────────┘   └────────┬───────────┘  │
│  └────────────┘                                     │              │
└───────────────────────────────────────────────────┬──────────────┘
                                                    │ USB / WiFi
┌───────────────────────────────────────────────────┴──────────────┐
│                ESP32 (Rust firmware, no_std)                      │
│  ┌──────────────────┐   ┌──────────────┐   ┌──────────────────┐  │
│  │  Recepción       │   │  Control PID │   │  Telemetría      │  │
│  │  de trayectoria  │   │  2 motores DC│   │  100 Hz          │  │
│  │  + Calibración   │   │  1 kHz       │   │  + NVS persist.  │  │
│  └──────────────────┘   └──────────────┘   └──────────────────┘  │
└──────────────────────────────────────────────────────────────────┘
```

## Requisitos

### Sistema operativo
- Linux (Ubuntu/Debian), Windows 10/11, o macOS

---

### 🐧 Linux (Ubuntu/Debian)

#### Dependencias del sistema
```bash
sudo apt-get install -y python3 python3-venv python3-full \
    libudev-dev pkg-config build-essential curl
```

#### Instalar Rust
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"
```

---

### 🪟 Windows 10/11

#### 1. Instalar Python

Descargar el instalador desde [python.org](https://www.python.org/downloads/) (versión 3.10 o superior).

> ⚠️ **Importante:** Durante la instalación, marcar la casilla **"Add Python to PATH"**.

Verificar la instalación:
```powershell
python --version
pip --version
```

#### 2. Instalar Rust

Descargar y ejecutar el instalador desde [rustup.rs](https://rustup.rs/).

O desde PowerShell:
```powershell
winget install Rustlang.Rustup
```

Verificar:
```powershell
rustc --version
cargo --version
```

#### 3. Opcional: Visual Studio Build Tools

Si al compilar `maturin develop` aparecen errores de linker, instalar
[Microsoft Visual C++ Build Tools](https://visualstudio.microsoft.com/visual-cpp-build-tools/).

En el instalador, seleccionar:
- **"Desarrollo para el escritorio con C++"** (Desktop development with C++)
- Asegurar que incluya **Windows SDK** y **MSVC v143**

> **Nota:** En la mayoría de los casos, Rust ya incluye su propio linker y no
> necesitas instalar Visual Studio. Solo es necesario si `maturin` falla con
> errores tipo `link.exe not found`.

---

## Instalación y ejecución

### 1. Clonar el repositorio

```bash
git clone <url-del-repo>
cd proyecto_academico
```

### 2. Crear entorno virtual e instalar dependencias Python

**Linux / macOS:**
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r backend/requirements.txt
```

**Windows (PowerShell / CMD):**
```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r backend\requirements.txt
```

### 3. Compilar e instalar el módulo Rust (opcional pero recomendado)

```bash
cd backend/rust_comm

# Compilar como extensión nativa de Python
maturin develop --release

cd ../..
```

> **Nota:** Si no compilas el módulo Rust, la aplicación usará automáticamente
> `pyserial` + `socket` como fallback. La funcionalidad es idéntica, solo
> cambia el rendimiento de la comunicación serie.

### 4. Ejecutar la aplicación

**Linux / macOS:**
```bash
.venv/bin/python backend/main.py
```

**Windows:**
```powershell
.venv\Scripts\python backend\main.py
```

## Uso de la interfaz

| Pestaña | Función |
|---------|---------|
| **1. Conexión** | Busca la ESP32 en red WiFi o USB. Permite configurar WiFi en la ESP32 vía USB. |
| **2. Calibración** | Asistente paso a paso para calibrar posición home, cero, límites angulares, torque de rozamiento y caracterización de curvas motor. |
| **3. Configuración** | Sliders interactivos con preview en tiempo real del trébol, torques y ángulos. |
| **4. Monitoreo** | Vista de la trayectoria de referencia (azul) y la trayectoria real (roja) recibida por telemetría. |

### Parámetros del trébol

| Parámetro | Rango | Descripción |
|-----------|-------|-------------|
| Pétalos | 3 – 7 | Número de hojas del trébol |
| Escala | 0.75 – 1.25 | Factor de escala de la figura |
| Velocidad | 1 – 10 cm/s | Velocidad tangencial de dibujado |
| Rotación | ±45° | Ángulo de rotación del trébol |
| Ciclos | 1 – 10 | Veces que se repite la trayectoria |

## Módulo de Calibración

El asistente de calibración (`backend/calibration/`) permite preparar el robot para
funcionar correctamente sin necesidad de conocer los detalles internos.

### Flujo de calibración

```
1. Verificar calibración previa (ESP32 NVS)
2. PWM rozamiento base → mover muy despacio hasta detectar movimiento
3. Límites angulares base → mover con PWM mínimo hasta stall mecánico
4. Cero base → el usuario indica cuál es la posición home (θ₁ = 175°)
5. PWM rozamiento codo → igual que base pero ya en posición home
6. Límites angulares codo
7. Cero codo (θ₂ = -128°)
8. Caracterización de curvas motor → barrido PWM 5%-100%
9. Ajuste función de transferencia G(s) = K/(τs + 1)
10. Guardar en memoria persistente (ESP32 NVS + archivo local)
```

### Tipos de calibración

| Tipo | Cuándo | Qué hace |
|------|--------|----------|
| **Completa** | Primera vez | Posición + caracterización de motores (FT) |
| **Solo posición** | Recalibración | Home, cero y límites (no toca curva motor) |

### Parámetros físicos

| Parámetro | Valor |
|-----------|-------|
| Eslabón 1 (base) | 195 mm centro a centro |
| Eslabón 2 (codo) | 260 mm centro a centro |
| Home θ₁ (base) | 190° (robot replegado a la izquierda) |
| Home θ₂ (codo) | -130° (plegado debajo de la línea media del trébol) |
| Motor base | 12V, 218:1, 30 kgf·cm |
| Motor codo | 24V, 270:1, 19 kgf·cm |

### Persistencia

- Los offsets de posición se guardan en **NVS de la ESP32** (sobrevive a apagados).
- Si la ESP32 no recuerda la posición al iniciar, se fuerza una recalibración.
- La caracterización de motores se realiza **una sola vez**.
- Backup local en `calibration_data.json`.

## Estructura del proyecto

```
proyecto_academico/
├── backend/
│   ├── main.py                      # Punto de entrada
│   ├── requirements.txt             # Dependencias Python
│   │
│   ├── gui/                         # Interfaz gráfica (customtkinter)
│   │   ├── main_window.py           # Ventana principal con tabs
│   │   ├── connection_dialog.py     # Flujo conexión red → USB
│   │   ├── calibration_wizard.py    # Asistente de calibración paso a paso
│   │   ├── config_panel.py          # Sliders + preview interactivo
│   │   └── trajectory_view.py       # Monitoreo en tiempo real
│   │
│   ├── calibration/                 # Calibración y setup inicial
│   │   ├── __init__.py              # Paquete de calibración
│   │   ├── calibration_manager.py   # Orquestador del flujo completo
│   │   ├── motor_identification.py  # Rozamiento, límites, curvas PWM
│   │   └── homing.py                # Home/cero de articulaciones
│   │
│   ├── trajectory/                  # Cálculo de trayectorias
│   │   ├── trefoil.py               # Generación del trébol + Bézier
│   │   ├── builder.py               # Orquestador de trayectoria
│   │   └── kinematics.py            # Cinemática inversa + dinámica
│   │
│   ├── communication/               # Protocolo PC ↔ ESP32
│   │   ├── __init__.py              # Trama binaria + CRC32
│   │   └── manager.py               # Gestor dual (Rust / pyserial)
│   │
│   ├── rust_comm/                   # Módulo Rust (PyO3) — PC
│   │   ├── Cargo.toml
│   │   ├── README.md
│   │   └── src/
│   │       ├── lib.rs               # Esp32Connection (Python API)
│   │       ├── protocol.rs          # Protocolo binario con CRC32
│   │       ├── discover.rs          # Auto-detección USB + UDP
│   │       └── serial_link.rs       # Enlace serie con buffering
│   │
│   └── esp32_firmware/              # Firmware ESP32 (Rust no_std)
│       ├── Cargo.toml
│       ├── README.md
│       └── src/
│           └── main.rs              # Control PID + comunicaciones
│
├── caracterizacion_trayectoria_trebol.m   # MATLAB original
├── modelo_motor_dc.m
├── caso_critico.m
├── datasheets/
└── esp32/                           # Código ESP32 (PlatformIO/C++)
```

## Licencia

Proyecto académico — Universidad, Servomecanismos.
