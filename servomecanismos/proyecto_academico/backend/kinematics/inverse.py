import numpy as np

# Dimensiones físicas del robot provenientes del modelo de motores (en metros)
L1 = 0.195  # 195 mm
L2 = 0.260  # 260 mm

def inverse_kinematics(x: float, y: float, elbow_up: bool = True):
    """
    Calcula la cinemática inversa de un robot planar 2R.
    
    Retorna:
        theta1, theta2 en grados para mover el efector final a (x, y) en metros.
    """
    # Verificación física simple del espacio de trabajo
    d_sq = x**2 + y**2
    d = np.sqrt(d_sq)
    
    if d > (L1 + L2) or d < abs(L1 - L2):
        # Punto fuera del espacio de trabajo, retorno por defecto seguro
        return 190.0, -130.0 # Posición de home
        
    # Ley de cosenos para theta2
    cos_theta2 = (d_sq - L1**2 - L2**2) / (2 * L1 * L2)
    cos_theta2 = np.clip(cos_theta2, -1.0, 1.0)
    
    if elbow_up:
        theta2 = -np.arccos(cos_theta2)
    else:
        theta2 = np.arccos(cos_theta2)
        
    # Ángulo theta1
    alpha = np.atan2(y, x)
    beta = np.atan2(L2 * np.sin(theta2), L1 + L2 * np.cos(theta2))
    theta1 = alpha - beta
    
    return np.degrees(theta1), np.degrees(theta2)
