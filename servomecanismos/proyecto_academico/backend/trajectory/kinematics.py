"""
Cinemática inversa y dinámica inversa para robot 2R (2 grados de libertad).

Implementa el modelo de Lagrange-Euler para un manipulador planar de dos
eslabones con masas concentradas en los centros de gravedad.
"""
import numpy as np
from scipy.signal import savgol_filter


# ── Parámetros físicos del robot ──────────────────────
l1 = 0.195          # Longitud eslabón 1 [m] (195 mm centro a centro)
l2 = 0.260          # Longitud eslabón 2 [m] (260 mm centro a centro)
m_link1 = 0.700      # Masa eslabón 1 [kg]
m_link2 = 0.600      # Masa eslabón 2 [kg]
cg1 = l1 / 2         # Centro de gravedad eslabón 1 [m]
cg2 = l2 / 2         # Centro de gravedad eslabón 2 [m]
m_motor2 = 0.300     # Masa del motor del codo (montado en eslabón 1)
m_tip = 0.050        # Masa de la punta/efector
g = 9.81             # Gravedad [m/s²]

# Masas y centros de gravedad equivalentes
m1_t = m_link1 + m_motor2
lc1_t = (m_link1 * cg1 + m_motor2 * l1) / m1_t

m2_t = m_link2 + m_tip
lc2_t = (m_link2 * cg2 + m_tip * l2) / m2_t

# Momentos de inercia equivalentes
I1_t = (1/12) * m_link1 * l1**2 + m_link1 * (lc1_t - cg1)**2 + m_motor2 * (l1 - lc1_t)**2
I2_t = (1/12) * m_link2 * l2**2 + m_link2 * (lc2_t - cg2)**2 + m_tip * (l2 - lc2_t)**2


def inverse_kinematics(
    x: np.ndarray,
    y: np.ndarray,
    th1_home: float,
    th2_home: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Cinemática inversa para robot 2R (configuración codo arriba).

    Args:
        x, y: Coordenadas cartesianas de la trayectoria.
        th1_home, th2_home: Ángulos iniciales (home) para el primer punto.

    Returns:
        (theta1, theta2): Arrays de ángulos articulares [rad].
    """
    n = len(x)
    theta1 = np.zeros(n)
    theta2 = np.zeros(n)

    for i in range(n):
        if i == 0:
            theta1[i] = th1_home
            theta2[i] = th2_home
            continue

        xi, yi = float(x[i]), float(y[i])

        D = (xi**2 + yi**2 - l1**2 - l2**2) / (2 * l1 * l2)
        D = np.clip(D, -1.0, 1.0)

        th2 = np.arctan2(-np.sqrt(1 - D**2), D)  # codo arriba
        th1 = np.arctan2(yi, xi) - np.arctan2(l2 * np.sin(th2), l1 + l2 * np.cos(th2))

        theta1[i] = th1
        theta2[i] = th2

    # Desenvolver ángulos
    theta1 = np.unwrap(theta1)
    theta2 = np.unwrap(theta2)

    return theta1, theta2


def inverse_dynamics(
    theta1: np.ndarray,
    theta2: np.ndarray,
    dt: float,
    poly_order: int = 3,
    window_size: int = 5,
) -> dict:
    """
    Dinámica inversa: calcula torques articulares a partir de la trayectoria.

    Args:
        theta1, theta2: Posiciones articulares [rad].
        dt: Paso de tiempo [s].
        poly_order, window_size: Parámetros del filtro Savitzky-Golay.

    Returns:
        Dict con velocidades, aceleraciones, torques, y métricas.
    """
    # Velocidades y aceleraciones
    dtheta1 = np.gradient(theta1, dt)
    dtheta2 = np.gradient(theta2, dt)

    ddtheta1 = np.gradient(dtheta1, dt)
    ddtheta2 = np.gradient(dtheta2, dt)

    # Filtrado Savitzky-Golay
    theta1_f = savgol_filter(theta1, window_size, poly_order)
    theta2_f = savgol_filter(theta2, window_size, poly_order)
    dtheta1_f = savgol_filter(dtheta1, window_size, poly_order)
    dtheta2_f = savgol_filter(dtheta2, window_size, poly_order)
    ddtheta1_f = savgol_filter(ddtheta1, window_size, poly_order)
    ddtheta2_f = savgol_filter(ddtheta2, window_size, poly_order)

    n = len(theta1)
    tau1 = np.zeros(n)
    tau2 = np.zeros(n)
    tau1_inercial = np.zeros(n)
    tau2_inercial = np.zeros(n)
    tau1_no_inercial = np.zeros(n)
    tau2_no_inercial = np.zeros(n)

    for i in range(n):
        th1 = theta1_f[i]
        th2 = theta2_f[i]
        d1 = dtheta1_f[i]
        d2 = dtheta2_f[i]
        dd1 = ddtheta1_f[i]
        dd2 = ddtheta2_f[i]

        # Matriz de inercia
        M11 = m1_t * lc1_t**2 + m2_t * (l1**2 + lc2_t**2 + 2 * l1 * lc2_t * np.cos(th2)) + I1_t + I2_t
        M12 = m2_t * (lc2_t**2 + l1 * lc2_t * np.cos(th2)) + I2_t
        M22 = m2_t * lc2_t**2 + I2_t

        # Términos de Coriolis y centrífugos
        h = -m2_t * l1 * lc2_t * np.sin(th2)
        C1 = h * (2 * d1 * d2 + d2**2)
        C2 = -h * d1**2

        # Gravedad
        G1 = (m1_t * lc1_t + m2_t * l1) * g * np.cos(th1) + m2_t * lc2_t * g * np.cos(th1 + th2)
        G2 = m2_t * lc2_t * g * np.cos(th1 + th2)

        tau1_inercial[i] = M11 * dd1 + M12 * dd2
        tau2_inercial[i] = M12 * dd1 + M22 * dd2
        tau1_no_inercial[i] = C1 + G1
        tau2_no_inercial[i] = C2 + G2

        tau1[i] = tau1_inercial[i] + tau1_no_inercial[i]
        tau2[i] = tau2_inercial[i] + tau2_no_inercial[i]

    results = {
        'theta1': theta1_f,
        'theta2': theta2_f,
        'dtheta1': dtheta1_f,
        'dtheta2': dtheta2_f,
        'ddtheta1': ddtheta1_f,
        'ddtheta2': ddtheta2_f,
        'tau1': tau1,
        'tau2': tau2,
        'tau1_inercial': tau1_inercial,
        'tau2_inercial': tau2_inercial,
        'tau1_no_inercial': tau1_no_inercial,
        'tau2_no_inercial': tau2_no_inercial,
    }
    return results


def forward_kinematics(
    theta1: np.ndarray,
    theta2: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Cinemática directa para el efector final."""
    x = l1 * np.cos(theta1) + l2 * np.cos(theta1 + theta2)
    y = l1 * np.sin(theta1) + l2 * np.sin(theta1 + theta2)
    return x, y
