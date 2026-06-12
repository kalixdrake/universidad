"""
Generación de trayectoria de trébol (trefoil) con N pétalos.

Port directo de caracterizacion_trayectoria_trebol.m
"""
import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter


def generate_trefoil_raw(
    L: float = 0.2,
    n_petals: int = 7,
    scale: float = 1.25,
    beta_deg: float = 45.0,
    num_points: int = 1000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Genera la trayectoria base del trébol en crudo (sin parametrización por
    velocidad ni desplazamiento).

    Retorna (phi_raw, x_raw, y_raw, r_raw).
    """
    r0 = (L / 2) * 0.75
    A = (L / 2) * 0.25
    beta = np.deg2rad(beta_deg)

    phi_raw = np.linspace(0, 2 * np.pi, num_points)
    r_raw = r0 + A * np.cos(n_petals * phi_raw)
    x_raw = scale * r_raw * np.cos(phi_raw + beta)
    y_raw = scale * r_raw * np.sin(phi_raw + beta)

    return phi_raw, x_raw, y_raw, r_raw


def parametrize_constant_speed(
    x_raw: np.ndarray,
    y_raw: np.ndarray,
    phi_raw: np.ndarray,
    v_const: float = 0.1,
    dt: float = 0.01,
    n_petals: int = 7,
    scale: float = 1.25,
    beta_deg: float = 45.0,
    L: float = 0.2,
    x_center: float = 0.2,
    y_center: float = 0.2,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Reparametriza la trayectoria para velocidad tangencial constante.

    Retorna (t, x_traj, y_traj, phi_t).
    """
    r0 = (L / 2) * 0.75
    A = (L / 2) * 0.25
    beta = np.deg2rad(beta_deg)

    # Longitud de arco acumulada
    dx = np.diff(x_raw)
    dy = np.diff(y_raw)
    ds = np.sqrt(dx**2 + dy**2)
    s_acum = np.concatenate([[0], np.cumsum(ds)])
    T_total = s_acum[-1] / v_const

    t = np.arange(0, T_total + dt, dt)
    s_t = v_const * t

    phi_t = interp1d(s_acum, phi_raw, kind='pchip', bounds_error=False, fill_value='extrapolate')(s_t)

    r_t = r0 + A * np.cos(n_petals * phi_t)
    x_traj = scale * r_t * np.cos(phi_t + beta) + x_center
    y_traj = scale * r_t * np.sin(phi_t + beta) + y_center

    # Cerrar la trayectoria
    x_traj[-1] = x_traj[0]
    y_traj[-1] = y_traj[0]

    return t, x_traj, y_traj, phi_t


def find_tangential_entry(
    x_traj: np.ndarray,
    y_traj: np.ndarray,
    x_home: float,
    y_home: float,
) -> tuple[np.ndarray, np.ndarray, int]:
    """
    Encuentra el punto de entrada a la trayectoria donde el vector desde
    home es más paralelo a la tangente local.  Reordena la trayectoria
    para empezar en ese punto.

    Retorna (x_reordered, y_reordered, idx_tangente).
    """
    dx_traj = np.gradient(x_traj)
    dy_traj = np.gradient(y_traj)
    norm_traj = np.sqrt(dx_traj**2 + dy_traj**2)
    tx = dx_traj / norm_traj
    ty = dy_traj / norm_traj

    hx = x_traj - x_home
    hy = y_traj - y_home
    norm_h = np.sqrt(hx**2 + hy**2)
    hx = hx / norm_h
    hy = hy / norm_h

    alineacion = hx * tx + hy * ty
    idx_tangente = int(np.argmax(alineacion))

    # Reordenar (sin duplicar el último punto)
    x_reordered = np.concatenate([x_traj[idx_tangente:-1], x_traj[:idx_tangente], [x_traj[idx_tangente]]])
    y_reordered = np.concatenate([y_traj[idx_tangente:-1], y_traj[:idx_tangente], [y_traj[idx_tangente]]])

    return x_reordered, y_reordered, idx_tangente


def generate_bezier_approach(
    x_start: float,
    y_start: float,
    x_target: float,
    y_target: float,
    x_target_next: float,
    y_target_next: float,
    v_const: float = 0.1,
    dt: float = 0.01,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Genera una curva de Bezier cúbica de aproximación desde (x_start, y_start)
    hasta (x_target, y_target) con entrada tangente a la trayectoria.

    Retorna (x_aprox, y_aprox).
    """
    dist_aprox = np.sqrt((x_target - x_start)**2 + (y_target - y_start)**2)

    P0x, P0y = x_start, y_start
    P3x, P3y = x_target, y_target

    # Vector tangente de entrada al trébol
    dx_llegada = x_target_next - x_target
    dy_llegada = y_target_next - y_target
    norm_llegada = np.sqrt(dx_llegada**2 + dy_llegada**2)
    tx_in = dx_llegada / norm_llegada
    ty_in = dy_llegada / norm_llegada

    L_ctrl = dist_aprox * 0.55

    # P2: antes del punto de llegada en dirección de la tangente
    P2x = P3x - L_ctrl * tx_in
    P2y = P3y - L_ctrl * ty_in

    # P1: arco visible hacia un lado
    vx = P3x - P0x
    vy = P3y - P0y
    norm_v = np.sqrt(vx**2 + vy**2)
    if norm_v > 1e-9:
        nx = -vy / norm_v
        ny = vx / norm_v
    else:
        nx, ny = 1.0, 0.0

    P1x = P0x + vx * 0.3 + nx * L_ctrl
    P1y = P0y + vy * 0.3 + ny * L_ctrl

    # Tiempo de aproximación
    T_aprox = (3 * L_ctrl) / v_const
    t_aprox = np.arange(0, T_aprox, dt)
    x_norm = t_aprox / T_aprox

    # Perfil suave: tau(0)=0, tau'(0)=0, tau''(0)=0, tau(1)=1, tau'(1)=0, tau''(1)=0
    tau = 3 * x_norm**5 - 8 * x_norm**4 + 6 * x_norm**3

    # Bezier cúbica
    x_aprox = (
        (1 - tau)**3 * P0x
        + 3 * (1 - tau)**2 * tau * P1x
        + 3 * (1 - tau) * tau**2 * P2x
        + tau**3 * P3x
    )
    y_aprox = (
        (1 - tau)**3 * P0y
        + 3 * (1 - tau)**2 * tau * P1y
        + 3 * (1 - tau) * tau**2 * P2y
        + tau**3 * P3y
    )

    # Quitar último punto (coincide con el primero del trébol)
    x_aprox = x_aprox[:-1]
    y_aprox = y_aprox[:-1]

    return x_aprox, y_aprox


def filter_savgol(
    data: np.ndarray,
    poly_order: int = 3,
    window_size: int = 5,
) -> np.ndarray:
    """Aplica filtro Savitzky-Golay."""
    return savgol_filter(data, window_size, poly_order)
