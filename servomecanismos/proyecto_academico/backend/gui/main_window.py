"""
Ventana principal de la aplicación Robot 2R.
Gestiona el flujo completo: conexión → configuración → cálculo → ejecución → monitoreo.
"""
import customtkinter as ctk
import logging

from .connection_dialog import ConnectionDialog
from .config_panel import ConfigPanel
from .trajectory_view import TrajectoryView

logger = logging.getLogger(__name__)

# Tema oscuro por defecto de customtkinter
ctk.set_appearance_mode("dark")


class MainWindow(ctk.CTk):
    """Ventana principal con navegación por tabs."""

    def __init__(self):
        super().__init__()

        self.title("Robot 2R — Control de Trayectoria")
        self.geometry("1200x800")
        self.minsize(1000, 650)

        # Estado
        self._trajectory_data: dict | None = None
        self._comm_mgr = None  # Se inicializa en ConnectionDialog

        # ── Layout ────────────────────────────────────
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # Tabview
        self.tabview = ctk.CTkTabview(self)
        self.tabview.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)

        self.tab_connect = self.tabview.add("1. Conexión")
        self.tab_config = self.tabview.add("2. Configuración")
        self.tab_monitor = self.tabview.add("3. Monitoreo")

        # Deshabilitar tabs hasta completar pasos
        self.tabview.set("1. Conexión")

        # ── Tab 1: Conexión ──────────────────────────
        self.connection_dialog = ConnectionDialog(
            self.tab_connect,
            on_connected=self._on_connected,
        )
        self.connection_dialog.pack(fill="both", expand=True, padx=5, pady=5)

        # ── Tab 2: Configuración (siempre visible) ────
        self.config_panel = ConfigPanel(
            self.tab_config,
            on_calculate=self._on_calculate,
        )
        self.config_panel.pack(fill="both", expand=True, padx=5, pady=5)

        # ── Tab 3: Monitoreo ──────────────────────────
        self.trajectory_view = TrajectoryView(
            self.tab_monitor,
            on_send=self._on_send_trajectory,
            on_stop=self._on_emergency_stop,
        )
        self.trajectory_view.pack(fill="both", expand=True, padx=5, pady=5)

        # ── Barra de estado ───────────────────────────
        self.status_var = ctk.StringVar(value="Listo — Conecte la ESP32 para comenzar")
        self.status_bar = ctk.CTkLabel(
            self, textvariable=self.status_var, anchor="w",
            fg_color="transparent", height=25,
        )
        self.status_bar.grid(row=1, column=0, sticky="ew", padx=10, pady=(0, 5))

        self.protocol("WM_DELETE_WINDOW", self._on_close)

    # ── Callbacks ─────────────────────────────────────
    def _on_connected(self, comm_mgr):
        """Callback llamado cuando se establece conexión con la ESP32."""
        self._comm_mgr = comm_mgr
        info = comm_mgr.info
        mode = "USB" if info.state.value == "usb" else "Red WiFi"
        self.status_var.set(f"Conectado por {mode} a ESP32 | Puerto: {info.port or info.address}")
        self.trajectory_view.set_connected(True)
        self.tabview.set("2. Configuración")

    def _on_calculate(self, trajectory: dict):
        """Callback cuando se presiona 'Calcular' en el panel de configuración."""
        self._trajectory_data = trajectory
        self.status_var.set(
            f"Trayectoria calculada: {len(trajectory['x_total'])} puntos, "
            f"dt={trajectory.get('dt', 0.01):.3f}s"
        )
        # Pasar datos a la vista de monitoreo
        self.trajectory_view.set_reference_trajectory(trajectory)
        self.tabview.set("3. Monitoreo")
        self._update_tab_state()

    def _on_send_trajectory(self):
        """Envía la trayectoria calculada a la ESP32."""
        if not self._comm_mgr:
            self.status_var.set("⚠ No hay ESP32 conectada. Vaya a la pestaña '1. Conexión' para conectar.")
            return
        if not self._trajectory_data:
            self.status_var.set("Error: No hay trayectoria calculada")
            return

        from ..communication import TrajectoryData

        dyn = self._trajectory_data.get('dynamics', {})
        cycles = self._trajectory_data.get('cycles', 1)
        dt = float(self._trajectory_data.get('dt', 0.01))

        theta1 = dyn.get('theta1', [])
        theta2 = dyn.get('theta2', [])
        omega1 = dyn.get('dtheta1', [])
        omega2 = dyn.get('dtheta2', [])

        traj = TrajectoryData(
            theta1=theta1.tolist() if hasattr(theta1, 'tolist') else list(theta1),
            theta2=theta2.tolist() if hasattr(theta2, 'tolist') else list(theta2),
            omega1=omega1.tolist() if hasattr(omega1, 'tolist') else list(omega1),
            omega2=omega2.tolist() if hasattr(omega2, 'tolist') else list(omega2),
            dt=dt,
            cycles=cycles,
        )

        self.status_var.set("Enviando trayectoria a la ESP32...")
        success = self._comm_mgr.send_trajectory(traj)
        if success:
            self._comm_mgr.send_start()
            self.status_var.set("Trayectoria enviada. ESP32 ejecutando...")
            # Configurar callback de telemetría
            self._comm_mgr.set_telemetry_callback(self.trajectory_view.on_telemetry)
        else:
            self.status_var.set("Error al enviar la trayectoria")

    def _on_emergency_stop(self):
        """Envía parada de emergencia a la ESP32."""
        if self._comm_mgr:
            self._comm_mgr.send_emergency_stop()
            self.status_var.set("⚠ PARADA DE EMERGENCIA enviada")

    def _update_tab_state(self):
        """Actualiza qué tabs están habilitados."""
        # TODO: customtkinter no soporta deshabilitar tabs nativamente,
        # se maneja por flujo de navegación.

    def _on_close(self):
        """Limpieza al cerrar la aplicación."""
        if self._comm_mgr:
            self._comm_mgr.close()
        self.destroy()
