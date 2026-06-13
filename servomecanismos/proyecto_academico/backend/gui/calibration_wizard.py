"""
Asistente gráfico de calibración del robot 2R.

Provee una interfaz paso a paso para:
1. Verificar estado previo de calibración
2. Calibrar posición home y cero de ambos motores
3. Caracterizar curvas velocidad/aceleración vs PWM
4. Guardar parámetros en memoria persistente

Diseñado para ser usable por cualquier persona, no solo el desarrollador.
"""
import customtkinter as ctk
import threading
import time
import logging
from typing import Optional, Callable

from ..calibration import (
    CalibrationManager,
    CalibrationState,
    CalibrationStep,
    CalibrationData,
    MotorParams,
)

logger = logging.getLogger(__name__)


class CalibrationWizard(ctk.CTkFrame):
    """
    Asistente de calibración completo con interfaz paso a paso.

    Uso:
        wizard = CalibrationWizard(parent, comm_manager=comm_mgr)
        wizard.pack(fill="both", expand=True)
    """

    def __init__(
        self,
        parent,
        comm_manager=None,
        on_calibration_complete: Optional[Callable[[CalibrationData], None]] = None,
    ):
        super().__init__(parent)
        self._comm = comm_manager
        self._on_complete = on_calibration_complete

        self._cal_mgr = CalibrationManager(comm_manager)
        self._cal_data: Optional[CalibrationData] = None
        self._running = False
        self._telemetry_cb: Optional[Callable[[], dict]] = None

        # ── Layout ────────────────────────────────────
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        # Título
        self.title_label = ctk.CTkLabel(
            self,
            text="🔧 Asistente de Calibración — Robot 2R",
            font=ctk.CTkFont(size=20, weight="bold"),
        )
        self.title_label.grid(row=0, column=0, padx=20, pady=(15, 5), sticky="w")

        # Panel principal (scrollable)
        self.main_frame = ctk.CTkScrollableFrame(self)
        self.main_frame.grid(row=1, column=0, padx=20, pady=10, sticky="nsew")
        self.main_frame.grid_columnconfigure(0, weight=1)

        # ── Sección de bienvenida ─────────────────────
        self._build_welcome_section()

        # ── Sección de estado ─────────────────────────
        self.status_section = ctk.CTkFrame(self.main_frame)
        self.status_section.grid(row=10, column=0, padx=10, pady=10, sticky="ew")
        self.status_section.grid_columnconfigure(0, weight=1)

        self.status_label = ctk.CTkLabel(
            self.status_section,
            text="Estado: Listo para comenzar",
            font=ctk.CTkFont(size=14),
            anchor="w",
        )
        self.status_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")

        # Barra de progreso
        self.progress_bar = ctk.CTkProgressBar(self.status_section)
        self.progress_bar.grid(row=1, column=0, padx=10, pady=5, sticky="ew")
        self.progress_bar.set(0)

        # Log de calibración
        self.log_text = ctk.CTkTextbox(self.status_section, height=120, font=ctk.CTkFont(size=12))
        self.log_text.grid(row=2, column=0, padx=10, pady=(10, 5), sticky="ew")

        # ── Botones de acción ─────────────────────────
        self.button_frame = ctk.CTkFrame(self.main_frame)
        self.button_frame.grid(row=20, column=0, padx=10, pady=15, sticky="ew")

        self.btn_check = ctk.CTkButton(
            self.button_frame,
            text="🔍 Verificar Calibración Existente",
            command=self._on_check_calibration,
            height=40,
            font=ctk.CTkFont(size=14),
        )
        self.btn_check.grid(row=0, column=0, padx=5, pady=5)

        self.btn_full = ctk.CTkButton(
            self.button_frame,
            text="⚙️ Calibración Completa (Posición + Motores)",
            command=self._on_full_calibration,
            height=40,
            fg_color="#2E7D32",
            hover_color="#1B5E20",
            font=ctk.CTkFont(size=14),
        )
        self.btn_full.grid(row=0, column=1, padx=5, pady=5)

        self.btn_position = ctk.CTkButton(
            self.button_frame,
            text="📍 Solo Calibración de Posición",
            command=self._on_position_calibration,
            height=40,
            fg_color="#1565C0",
            hover_color="#0D47A1",
            font=ctk.CTkFont(size=14),
        )
        self.btn_position.grid(row=0, column=2, padx=5, pady=5)

        self.btn_stop = ctk.CTkButton(
            self.button_frame,
            text="🛑 PARADA DE EMERGENCIA",
            command=self._on_emergency_stop,
            height=40,
            fg_color="#C62828",
            hover_color="#8E0000",
            font=ctk.CTkFont(size=14),
        )
        self.btn_stop.grid(row=0, column=3, padx=5, pady=5)

        # ── Opciones avanzadas ────────────────────────
        self.advanced_frame = ctk.CTkFrame(self.main_frame)
        self.advanced_frame.grid(row=15, column=0, padx=10, pady=10, sticky="ew")
        self.advanced_frame.grid_columnconfigure((0, 1), weight=1)

        ctk.CTkLabel(
            self.advanced_frame,
            text="Ángulo Home Personalizado (opcional):",
            font=ctk.CTkFont(size=13, weight="bold"),
        ).grid(row=0, column=0, columnspan=2, padx=10, pady=(5, 0), sticky="w")

        ctk.CTkLabel(self.advanced_frame, text="Base θ₁ (°):").grid(
            row=1, column=0, padx=10, pady=2, sticky="e"
        )
        self.home_th1_var = ctk.StringVar(value="190.0")
        self.home_th1_entry = ctk.CTkEntry(self.advanced_frame, width=100, textvariable=self.home_th1_var)
        self.home_th1_entry.grid(row=1, column=1, padx=10, pady=2, sticky="w")

        ctk.CTkLabel(self.advanced_frame, text="Codo θ₂ (°):").grid(
            row=2, column=0, padx=10, pady=2, sticky="e"
        )
        self.home_th2_var = ctk.StringVar(value="-130.0")
        self.home_th2_entry = ctk.CTkEntry(self.advanced_frame, width=100, textvariable=self.home_th2_var)
        self.home_th2_entry.grid(row=2, column=1, padx=10, pady=2, sticky="w")

        # ── Info de parámetros ────────────────────────
        self.params_frame = ctk.CTkFrame(self.main_frame)
        self.params_frame.grid(row=16, column=0, padx=10, pady=10, sticky="ew")
        self.params_frame.grid_columnconfigure((0, 1), weight=1)

        self.params_text = ctk.CTkTextbox(
            self.params_frame, height=180,
            font=ctk.CTkFont(size=12, family="monospace"),
        )
        self.params_text.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.params_text.insert("1.0", "Parámetros de calibración no disponibles.\nEjecute una calibración primero.")

        # ── Inicializar callbacks del manager ──────────
        self._cal_mgr.set_step_callback(self._on_calibration_step)
        self._cal_mgr.set_state_callback(self._on_calibration_state)

    # ── Construcción de UI ────────────────────────────

    def _build_welcome_section(self):
        """Sección de bienvenida con instrucciones."""
        welcome = ctk.CTkFrame(self.main_frame)
        welcome.grid(row=0, column=0, padx=10, pady=10, sticky="ew")
        welcome.grid_columnconfigure(0, weight=1)

        instructions = (
            "🎯 Bienvenido al Asistente de Calibración\n\n"
            "Este proceso preparará el robot 2R para funcionar correctamente.\n\n"
            "📋 Flujo de calibración:\n"
            "  ┌─────────────────────────────────────────────────┐\n"
            "  │ 1. Verificar calibración previa (ESP32 NVS)     │\n"
            "  │ 2. Calibrar posición Home y cero de cada motor  │\n"
            "  │ 3. Encontrar torque de rozamiento (PWM mínimo)  │\n"
            "  │ 4. Buscar límites angulares de cada articulación│\n"
            "  │ 5. Caracterizar curva velocidad/PWM de motores  │\n"
            "  │ 6. Ajustar función de transferencia G(s)=K/(τs+1)│\n"
            "  │ 7. Guardar todo en memoria persistente          │\n"
            "  └─────────────────────────────────────────────────┘\n\n"
            "⚠️  Recomendaciones:\n"
            "  • Asegúrese de que el robot tenga espacio para moverse.\n"
            "  • La calibración de motores solo se necesita UNA vez.\n"
            "  • La calibración de posición se puede repetir si es necesario.\n"
            "  • Use PARADA DE EMERGENCIA si algo sale mal.\n\n"
            f"📐 Home por defecto: Base θ₁ = {self._cal_mgr.home_th1_deg}°, "
            f"Codo θ₂ = {self._cal_mgr.home_th2_deg}°\n"
            f"   (Robot replegado, debajo del trébol y a la izquierda)"
        )

        self.welcome_label = ctk.CTkLabel(
            welcome,
            text=instructions,
            font=ctk.CTkFont(size=13),
            justify="left",
            anchor="w",
        )
        self.welcome_label.grid(row=0, column=0, padx=15, pady=15, sticky="w")

    # ── Setters ───────────────────────────────────────

    def set_comm_manager(self, comm_manager):
        """Actualiza el gestor de comunicación."""
        self._comm = comm_manager
        self._cal_mgr = CalibrationManager(comm_manager)
        self._cal_mgr.set_step_callback(self._on_calibration_step)
        self._cal_mgr.set_state_callback(self._on_calibration_state)

    def set_telemetry_callback(self, callback: Callable[[], dict]):
        """Establece callback para obtener telemetría en tiempo real."""
        self._telemetry_cb = callback

    # ── Acciones de botones ───────────────────────────

    def _on_check_calibration(self):
        """Verifica si existe calibración previa."""
        self._log("🔍 Verificando calibración existente...")
        state = self._cal_mgr.check_state()

        if state == CalibrationState.FULLY_CALIBRATED:
            self._log("✅ Robot COMPLETAMENTE calibrado.")
            self._log("   - Posición calibrada: ✓")
            self._log("   - Motores caracterizados: ✓")
            self._show_params()
        elif state == CalibrationState.POSITION_ONLY:
            self._log("⚠️  Solo posición calibrada. Falta caracterización de motores.")
            self._show_params()
        elif state == CalibrationState.NEEDS_RECALIBRATION:
            self._log("⚠️  Datos de calibración corruptos. Se necesita recalibración.")
        else:
            self._log("❌ No hay calibración previa. Ejecute calibración completa.")

    def _on_full_calibration(self):
        """Ejecuta calibración completa en un hilo separado."""
        if self._running:
            self._log("⚠️  Ya hay una calibración en curso.")
            return

        self._running = True
        self._log("⚙️ Iniciando calibración COMPLETA...")
        self._log("   Esto incluye: posición, límites y caracterización de motores.")
        self._set_buttons_state("disabled")

        th1 = float(self.home_th1_var.get() or "190.0")
        th2 = float(self.home_th2_var.get() or "-130.0")

        def run():
            try:
                self._cal_data = self._cal_mgr.run_full_calibration(
                    telemetry_callback=self._telemetry_cb,
                    user_home_th1_deg=th1,
                    user_home_th2_deg=th2,
                )
                self.after(0, lambda: self._on_calibration_finished(success=True))
            except Exception as e:
                logger.error(f"Error en calibración: {e}")
                self.after(0, lambda: self._log(f"❌ Error: {e}"))
                self.after(0, lambda: self._on_calibration_finished(success=False))

        threading.Thread(target=run, daemon=True).start()

    def _on_position_calibration(self):
        """Ejecuta solo calibración de posición."""
        if self._running:
            self._log("⚠️  Ya hay una calibración en curso.")
            return

        self._running = True
        self._log("📍 Iniciando calibración de POSICIÓN...")
        self._log("   Esto incluye: home, cero y límites de ambos motores.")
        self._set_buttons_state("disabled")

        th1 = float(self.home_th1_var.get() or "190.0")
        th2 = float(self.home_th2_var.get() or "-130.0")

        def run():
            try:
                self._cal_data = self._cal_mgr.run_position_only_calibration(
                    telemetry_callback=self._telemetry_cb,
                    user_home_th1_deg=th1,
                    user_home_th2_deg=th2,
                )
                self.after(0, lambda: self._on_calibration_finished(success=True))
            except Exception as e:
                logger.error(f"Error en calibración de posición: {e}")
                self.after(0, lambda: self._log(f"❌ Error: {e}"))
                self.after(0, lambda: self._on_calibration_finished(success=False))

        threading.Thread(target=run, daemon=True).start()

    def _on_emergency_stop(self):
        """Parada de emergencia."""
        self._log("🛑 ¡PARADA DE EMERGENCIA!")
        self._running = False
        if self._comm:
            self._comm.stop_all_motors()
            self._comm.send_emergency_stop()
        self._set_buttons_state("normal")

    # ── Callbacks del manager ─────────────────────────

    def _on_calibration_step(self, step: CalibrationStep, message: str, progress: float):
        """Recibe actualizaciones de progreso del CalibrationManager."""
        self.after(0, lambda: self._update_progress(step, message, progress))

    def _on_calibration_state(self, state: CalibrationState):
        """Recibe cambios de estado."""
        state_names = {
            CalibrationState.NOT_CALIBRATED: "No calibrado",
            CalibrationState.NEEDS_RECALIBRATION: "Necesita recalibración",
            CalibrationState.POSITION_ONLY: "Solo posición calibrada",
            CalibrationState.FULLY_CALIBRATED: "Completamente calibrado",
            CalibrationState.IN_PROGRESS: "Calibración en curso...",
            CalibrationState.ERROR: "Error",
        }
        name = state_names.get(state, str(state))
        self.after(0, lambda: self.status_label.configure(text=f"Estado: {name}"))

    def _update_progress(self, step: CalibrationStep, message: str, progress: float):
        """Actualiza barra de progreso y log."""
        self.progress_bar.set(min(progress, 1.0))
        self._log(f"  [{step.name}] {message}")

    # ── Finalización ──────────────────────────────────

    def _on_calibration_finished(self, success: bool):
        """Maneja la finalización de la calibración."""
        self._running = False
        self._set_buttons_state("normal")

        if success and self._cal_data:
            self._log("✅ ¡Calibración completada exitosamente!")
            self._show_params()

            if self._cal_data.is_motor_characterized:
                self._log("   📊 Función de transferencia:")
                self._log(f"      Motor base: K={self._cal_data.motor1.K:.3f}, τ={self._cal_data.motor1.tau:.3f}s")
                self._log(f"      Motor codo: K={self._cal_data.motor2.K:.3f}, τ={self._cal_data.motor2.tau:.3f}s")

            if self._on_complete:
                self._on_complete(self._cal_data)
        else:
            self._log("❌ La calibración NO se completó correctamente.")

    def _show_params(self):
        """Muestra los parámetros de calibración en el panel de info."""
        if not self._cal_data:
            return

        data = self._cal_data
        m1 = data.motor1
        m2 = data.motor2

        text = (
            "══════════════════════════════════════════\n"
            "       PARÁMETROS DE CALIBRACIÓN          \n"
            "══════════════════════════════════════════\n\n"
            "📐 POSICIÓN HOME:\n"
            f"   θ₁ home offset:  {data.home_th1_offset*180/3.14159:.2f}°\n"
            f"   θ₂ home offset:  {data.home_th2_offset*180/3.14159:.2f}°\n\n"
            "⚙️  MOTOR BASE (M1):\n"
            f"   PWM rozamiento:  {m1.friction_pwm*100:.1f}%\n"
            f"   θ₁ min:          {m1.angle_min*180/3.14159:.1f}°\n"
            f"   θ₁ max:          {m1.angle_max*180/3.14159:.1f}°\n"
            f"   FT: G(s) = {m1.K:.3f} / ({m1.tau:.3f}s + 1)\n"
            f"   R² ajuste:       {m1.r_squared:.4f}\n"
            f"   Puntos calib:    {len(m1.pwm_points)}\n\n"
            "⚙️  MOTOR CODO (M2):\n"
            f"   PWM rozamiento:  {m2.friction_pwm*100:.1f}%\n"
            f"   θ₂ min:          {m2.angle_min*180/3.14159:.1f}°\n"
            f"   θ₂ max:          {m2.angle_max*180/3.14159:.1f}°\n"
            f"   FT: G(s) = {m2.K:.3f} / ({m2.tau:.3f}s + 1)\n"
            f"   R² ajuste:       {m2.r_squared:.4f}\n"
            f"   Puntos calib:    {len(m2.pwm_points)}\n\n"
            "══════════════════════════════════════════\n"
        )

        if data.position_calibrated_at:
            text += f"\n📅 Posición calibrada:  {data.position_calibrated_at[:19]}"
        if data.motor_characterized_at:
            text += f"\n📅 Motores caracterizados: {data.motor_characterized_at[:19]}"

        self.params_text.delete("1.0", "end")
        self.params_text.insert("1.0", text)

    # ── Utilidades ────────────────────────────────────

    def _log(self, message: str):
        """Añade mensaje al log."""
        self.log_text.insert("end", message + "\n")
        self.log_text.see("end")

    def _set_buttons_state(self, state: str):
        """Habilita/deshabilita botones."""
        for btn in [self.btn_check, self.btn_full, self.btn_position, self.btn_stop]:
            try:
                btn.configure(state=state)
            except Exception:
                pass

    # ── Acceso a datos ────────────────────────────────

    def get_calibration_data(self) -> Optional[CalibrationData]:
        """Retorna los datos de calibración actuales."""
        return self._cal_data
