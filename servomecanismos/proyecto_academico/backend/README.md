# Requisitos del Sistema - Backend de Control PC

Este subproyecto contiene el software que se ejecuta en la computadora para controlar el Robot planar de 2 grados de libertad (2R), generar trayectorias, y comunicarse con la ESP32.

## Requisitos y Módulos del Sistema

Para ejecutar este software, el computador necesita:

1. **Python 3.10 o superior**
2. **Bibliotecas de Python**:
   - `customtkinter`: Para la interfaz gráfica de usuario (GUI) moderna.
   - `numpy`: Para manipulación de matrices y cálculo de trayectorias.
   - `scipy`: Para optimización y funciones matemáticas avanzadas.
   - `pyserial`: Para comunicarse por puerto serie (USB) con la ESP32.
   - `matplotlib`: Para previsualizar las trayectorias calculadas en la PC.

3. **Herramientas de Rust (en la PC)**:
   - **Maturin**: `pip install maturin` (herramienta de Python para compilar módulos de Rust como extensiones nativas binarias para Python).
   - **Rustup & Cargo**: La cadena de herramientas de Rust, requerida si deseas compilar la biblioteca nativa y optimizada de comunicación `rust_comm` para la PC.
   - **espflash**: `cargo install espflash` (herramienta CLI de Rust para interactuar, programar/flashear y monitorear dispositivos ESP32 por línea de comandos).

---

## Preparación del Entorno (Linux / macOS / Windows)

Se recomienda encarecidamente utilizar un entorno virtual de Python (`venv`) para no interferir con las dependencias del sistema global.

### 1. Activar el Entorno Virtual

Si ya existe un entorno virtual `.venv` en la raíz del proyecto académico:

**En Linux / macOS:**
```bash
source ../.venv/bin/activate
```

**En Windows (PowerShell):**
```powershell
..\.venv\Scripts\Activate.ps1
```

O si deseas crear uno nuevo en esta carpeta `backend`:
```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 2. Instalar los Requisitos

Instala todas las dependencias necesarias ejecutando:
```bash
pip install -r requirements.txt
```

---

## Estructura del Proyecto

- `main.py`: Punto de entrada del programa de control en PC.
- `trajectory/`: Generación matemática de trayectorias (ej. Trébol de N-pétalos, Curvas de Bézier).
- `kinematics/`: Cálculos de Cinemática Inversa y Directa para el robot planar 2R.
- `communication/`: Protocolo de comunicación (Serial/WiFi) con detección de errores (CRC32).
- `calibration/`: Rutinas de identificación del motor, autocalibración y homing.

---

## Cómo Ejecutar El Esqueleto

Una vez instalado el entorno virtual y las dependencias, ejecuta:
```bash
python main.py
```
