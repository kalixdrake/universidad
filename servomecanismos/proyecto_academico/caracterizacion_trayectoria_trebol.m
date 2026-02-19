%% Parámetros del Proyecto 
L = 0.20;           % Lado mínimo de 20cm
n = 7;              % Número de hojas (n=4 en la figura) 
v_const = 0.01;      % Velocidad entre 1 y 10 cm/s
S = 1.00;           % Escala inicial
beta = deg2rad(0);  % Rotación inicial (+/- 45°)

% Definición geométrica para centrar en el origen
r0 = (L/2) * 0.75;  
A = (L/2) * 0.25;   % r_max = r0 + A = L/2

%% 1. Generación de la Trayectoria Base (Garantizar Cierre)
% Usamos una resolución alta de phi para calcular la longitud de arco
phi_raw = linspace(0, 2*pi, 1000); 
r_raw = r0 + A*cos(n * phi_raw);
x_raw = r_raw .* cos(phi_raw);
y_raw = r_raw .* sin(phi_raw);

%% 2. Parametrización para Velocidad Constante
% Calculamos la distancia acumulada (arco)
dx = diff(x_raw); dy = diff(y_raw);
ds = sqrt(dx.^2 + dy.^2);
s_acumulada = [0, cumsum(ds)];
perimetro_total = s_acumulada(end);

% Definimos el tiempo total basado en la velocidad deseada
T_total = perimetro_total / v_const;
dt = 0.01; % Paso de tiempo de simulación
t = 0:dt:T_total;

% Interpolamos para obtener s(t) lineal (velocidad constante)
s_t = v_const * t; 

% Interpolamos phi para que corresponda a esas distancias s_t
phi_t = interp1(s_acumulada, phi_raw, s_t, 'pchip', 'extrap');

%% 3. Cálculo de Coordenadas Finales (Centradas en el Origen)
r_t = r0 + A*cos(n * phi_t);
x_traj = S * r_t .* cos(phi_t + beta);
y_traj = S * r_t .* sin(phi_t + beta);

% Asegurar cierre manual del último punto con el primero
x_traj(end) = x_traj(1);
y_traj(end) = y_traj(1);

%% 4. Fase de Aproximación 
% El robot inicia en un punto "recogido" a la izquierda [cite: 20]
x_inicio = -L/2; y_inicio = 0; 
x_aprox = linspace(x_inicio, x_traj(1), 50);
y_aprox = linspace(y_inicio, y_traj(1), 50);

%% Visualización
figure(1); clf; hold on;
plot(x_aprox, y_aprox, 'r--', 'DisplayName', 'Aproximación');
plot(x_traj, y_traj, 'b', 'LineWidth', 1.5, 'DisplayName', 'Trébol');
plot(0, 0, 'kx', 'MarkerSize', 10, 'DisplayName', 'Centro (Origen)');
axis equal; grid on; legend;
title(['Trayectoria de Trébol Estilizado (n=' num2str(n) ')']);
xlabel('X [m]'); ylabel('Y [m]');

%% Parámetros Físicos del Robot (en metros)
l1 = 0.15; % Longitud eslabón 1 [cite: 7]
l2 = 0.15; % Longitud eslabón 2 [cite: 6]

% Inicialización de vectores de ángulos
theta1 = zeros(size(x_traj));
theta2 = zeros(size(y_traj));

%% Cálculo de Cinemática Inversa (IK)
for i = 1:length(x_traj)
    px = x_traj(i);
    py = y_traj(i);
    
    % Cálculo de cos(theta2)
    cos_th2 = (px^2 + py^2 - l1^2 - l2^2) / (2 * l1 * l2);
    
    % Verificación de alcance (evitar errores numéricos fuera de [-1, 1])
    cos_th2 = max(min(cos_th2, 1), -1);
    
    % Elegimos "Codo arriba" (positivo) o "Codo abajo" (negativo)
    % Se recomienda mantener el mismo signo para toda la trayectoria
    theta2(i) = acos(cos_th2); 
    
    % Cálculo de theta1
    theta1(i) = atan2(py, px) - atan2(l2 * sin(theta2(i)), l1 + l2 * cos(theta2(i)));
end
%% 

%% Visualización de los Ángulos de los Motores
figure(2);
subplot(2,1,1);
plot(t, rad2deg(theta1), 'r', 'LineWidth', 1.5);
title('Ángulo Motor 1 (\theta_1)'); ylabel('Grados [°]'); grid on;

subplot(2,1,2);
plot(t, rad2deg(theta2), 'b', 'LineWidth', 1.5);
title('Ángulo Motor 2 (\theta_2)'); ylabel('Grados [°]'); xlabel('Tiempo [s]'); grid on;

%% Simulación del Movimiento del Robot
figure(3);
for i = 1:20:length(x_traj) % Dibujar cada 20 pasos para rapidez
    % Coordenadas de las articulaciones
    x0 = 0; y0 = 0;
    x1 = l1 * cos(theta1(i));
    y1 = l1 * sin(theta1(i));
    x2 = x_traj(i);
    y2 = y_traj(i);
    
    clf;
    hold on;
    plot(x_traj, y_traj, 'k:'); % Trayectoria deseada
    plot([x0 x1], [y0 y1], 'r-o', 'LineWidth', 2); % Eslabón 1
    plot([x1 x2], [y1 y2], 'b-o', 'LineWidth', 2); % Eslabón 2
    axis equal; grid on;
    axis([-0.3 0.3 -0.1 0.4]);
    title('Simulación Robot 2R siguiendo Trébol');
    pause(0.01);
end