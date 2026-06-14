import numpy as np

def generate_clover_trajectory(petals: int = 3, num_points: int = 200):
    """
    Genera una trayectoria en forma de trébol (pétalos) en coordenadas cartesianas (x, y).
    
    Parámetros:
        petals: Número de pétalos del trébol.
        num_points: Resolución de la trayectoria.
    Retorna:
        Matriz de puntos (x, y) de forma (num_points, 2).
    """
    # Esqueleto matemático simple (coordenadas polares r = cos(k * theta))
    theta = np.linspace(0, 2 * np.pi, num_points)
    r = 0.1 * np.cos(petals * theta) + 0.15 # Escalado y offset conveniente para el robot
    
    x = r * np.cos(theta)
    y = r * np.sin(theta) + 0.20 # Desplazamiento hacia el espacio de trabajo del robot
    
    return np.column_stack((x, y))
