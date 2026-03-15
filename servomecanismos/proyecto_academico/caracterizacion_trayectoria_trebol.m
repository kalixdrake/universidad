%% CARACTERIZACIÓN, CINEMÁTICA Y DINÁMICA DE ROBOT 2R (OPTIMIZADO PARA ESP32)
clear; clc; close all;

%% 1. Parámetros del Proyecto 
L = 0.20;           % Lado
n = 7;              % Número de hojas
v_const = 0.05;      % Velocidad deseada
S = 1.25;           % Escala
beta = deg2rad(30); % Rotación

% Longitud de eslabones, masas y centros de gravedad
l1 = 0.18; m_link1 = 0.700; cg1 = l1 / 2; % cg1 es el centro de gravedad
l2 = 0.18; m_link2 = 0.600; cg2 = l2 / 2; % Igual con cg2
m_motor2 = 0.300; m_tip = 0.050;
g = 9.81;

%% 2. Generación de la Trayectoria Base del Trébol
r0 = (L/2) * 0.75;  
A = (L/2) * 0.25;   

phi_raw = linspace(0, 2*pi, 1000); 
r_raw = r0 + A*cos(n * phi_raw);
x_raw = r_raw .* cos(phi_raw);
y_raw = r_raw .* sin(phi_raw);

% Parametrización para velocidad constante
dx = diff(x_raw); dy = diff(y_raw);
ds = sqrt(dx.^2 + dy.^2);
s_acumulada = [0, cumsum(ds)];
T_total = s_acumulada(end) / v_const;
dt = 0.01; % Paso de tiempo de simulación
t = 0:dt:T_total;

s_t = v_const * t; 
phi_t = interp1(s_acumulada, phi_raw, s_t, 'pchip', 'extrap');

r_t = r0 + A*cos(n * phi_t);
x_traj = S * r_t .* cos(phi_t + beta);
y_traj = S * r_t .* sin(phi_t + beta);
x_traj(end) = x_traj(1); y_traj(end) = y_traj(1);

% Desplazamiento del Origen al espacio de trabajo seguro
X_centro = 0.1;  Y_centro = 0.19; 
x_traj_desp = x_traj + X_centro;
y_traj_desp = y_traj + Y_centro;

%% 3. Posición "Home" (Codo Arriba)
disp('1. Calculando posición Home óptima...');
th1_home = pi/2; % Eslabón 1 vertical
Y_codo = l1 * sin(th1_home);

th2_test = linspace(-pi, 0, 180); % Solo codo arriba
Y_efector_test = Y_codo + l2 * sin(th1_home + th2_test);

posiciones_validas = (Y_efector_test >= 0.02) & (Y_efector_test <= Y_centro);

% Calcular torques estáticos aproximados para elegir el menor esfuerzo (
m2_t = m_link2 + m_tip;
lc2_t = (m_link2 * cg2 + m_tip * l2) / m2_t; % Usando el nuevo cg2
G1_test = m2_t * lc2_t * g * cos(th1_home + th2_test);
G2_test = m2_t * lc2_t * g * cos(th1_home + th2_test);
Torque_Estacionario = abs(G1_test) + abs(G2_test);
Torque_Estacionario(~posiciones_validas) = inf;

[~, idx] = min(Torque_Estacionario);
th2_home = th2_test(idx);

x_inicio = l1 * cos(th1_home) + l2 * cos(th1_home + th2_home);
y_inicio = l1 * sin(th1_home) + l2 * sin(th1_home + th2_home);

%% 4. Punto de inicio tangencial (NUEVO)
disp('2. Buscando punto de entrada tangencial...');
% 1. Encontrar el vector tangente en cada punto de la trayectoria
dx_traj = gradient(x_traj_desp);
dy_traj = gradient(y_traj_desp);
norm_traj = sqrt(dx_traj.^2 + dy_traj.^2);
tx = dx_traj ./ norm_traj; 
ty = dy_traj ./ norm_traj;

% 2. Encontrar el vector desde Home hacia cada punto
hx = x_traj_desp - x_inicio;
hy = y_traj_desp - y_inicio;
norm_h = sqrt(hx.^2 + hy.^2);
hx = hx ./ norm_h;
hy = hy ./ norm_h;

% 3. Producto punto para forzar que sea paralela a la tangente
alineacion = hx .* tx + hy .* ty;
[~, idx_tangente] = max(alineacion); % El valor más cercano a 1

x_temp = x_traj_desp(1:end-1); y_temp = y_traj_desp(1:end-1);
x_traj_opt = [x_temp(idx_tangente:end), x_temp(1:idx_tangente-1)];
y_traj_opt = [y_temp(idx_tangente:end), y_temp(1:idx_tangente-1)];

x_traj_desp = [x_traj_opt, x_traj_opt(1)];
y_traj_desp = [y_traj_opt, y_traj_opt(1)];

%% 5. Fase de Aproximación
disp('3. Generando aproximación ...');
dist_aprox = sqrt((x_traj_desp(1) - x_inicio)^2 + (y_traj_desp(1) - y_inicio)^2);

% Perfil para acelerar desde v=0 hasta v=v_const en la distancia dada
T_aprox = 2 * dist_aprox / v_const; 
t_aprox = 0:dt:T_aprox;

% Perfil cuadrático (para evitar saltos de torque abruptos)
tau = (t_aprox / T_aprox).^2; 

x_aprox = x_inicio + (x_traj_desp(1) - x_inicio) * tau;
y_aprox = y_inicio + (y_traj_desp(1) - y_inicio) * tau;
x_aprox(end) = []; y_aprox(end) = [];

x_total = [x_aprox, x_traj_desp];
y_total = [y_aprox, y_traj_desp];
num_puntos = length(x_total);

%% 6. Cinemática Inversa Vectorizada (Codo Arriba)
disp('4. Calculando Cinemática Inversa...');
theta1_raw = zeros(1, num_puntos);
theta2_raw = zeros(1, num_puntos);

for i = 1:num_puntos
    x = x_total(i); y = y_total(i);
    D = (x^2 + y^2 - l1^2 - l2^2) / (2 * l1 * l2);
    if D > 1, D = 1; elseif D < -1, D = -1; end
    
    th2 = atan2(-sqrt(1 - D^2), D);  
    th1 = atan2(y, x) - atan2(l2 * sin(th2), l1 + l2 * cos(th2));
    
    theta1_raw(i) = th1; theta2_raw(i) = th2;
end
theta1_raw = unwrap(theta1_raw);
theta2_raw = unwrap(theta2_raw);

%% 7. Filtrado SavGol (Básicamente no se usa ahora, pero será necesario para procesar las respuestas del microcontrolador)
disp('5. Optimizando perfiles para motores (Filtrado Savitzky-Golay)...');
orden_pol = 3;       
window_size = 5;    

theta1 = sgolayfilt(theta1_raw, orden_pol, window_size);
theta2 = sgolayfilt(theta2_raw, orden_pol, window_size);

dtheta1 = sgolayfilt(gradient(theta1, dt), orden_pol, window_size);
dtheta2 = sgolayfilt(gradient(theta2, dt), orden_pol, window_size);

ddtheta1 = sgolayfilt(gradient(dtheta1, dt), orden_pol, window_size);
ddtheta2 = sgolayfilt(gradient(dtheta2, dt), orden_pol, window_size);

t_total_vec = (0:num_puntos-1) * dt;

%% 8. Que tanto se aleja el filtrado de el ideal?
x_ideal = l1 * cos(theta1_raw) + l2 * cos(theta1_raw + theta2_raw);
y_ideal = l1 * sin(theta1_raw) + l2 * sin(theta1_raw + theta2_raw);

x_real = l1 * cos(theta1) + l2 * cos(theta1 + theta2);
y_real = l1 * sin(theta1) + l2 * sin(theta1 + theta2);

figure(1); clf; hold on; axis equal; grid on;
plot(x_ideal, y_ideal, 'Color', [0.8 0.8 0.8], 'LineWidth', 3, 'DisplayName', 'Trébol Matemático Ideal');
plot(x_real, y_real, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trébol Suavizado (Motores)');
title('Comparación: Trayectoria Ideal vs Suavizada');
legend('Location', 'best');

%% 9. Dinámica Inversa y Cálculo de Torques (ACTUALIZADA CON cg1 y cg2)
disp('6. Calculando Dinámica Inversa...');
m1_t = m_link1 + m_motor2;
lc1_t = (m_link1 * cg1 + m_motor2 * l1) / m1_t; 
% NOTA: Si tienes el Izz de CAD, reemplaza el término (1/12)*m*l^2 por ese Izz
I1_t = (1/12)*m_link1*l1^2 + m_link1*(lc1_t - cg1)^2 + m_motor2*(l1 - lc1_t)^2; 

I2_t = (1/12)*m_link2*l2^2 + m_link2*(lc2_t - cg2)^2 + m_tip*(l2 - lc2_t)^2;

tau1 = zeros(1, num_puntos); tau2 = zeros(1, num_puntos);

for i = 1:num_puntos
    th1 = theta1(i); th2 = theta2(i);
    d1 = dtheta1(i); d2 = dtheta2(i);
    dd1 = ddtheta1(i); dd2 = ddtheta2(i);
    
    M11 = m1_t*lc1_t^2 + m2_t*(l1^2 + lc2_t^2 + 2*l1*lc2_t*cos(th2)) + I1_t + I2_t;
    M12 = m2_t*(lc2_t^2 + l1*lc2_t*cos(th2)) + I2_t;
    M22 = m2_t*lc2_t^2 + I2_t;
    
    h_coriolis = -m2_t*l1*lc2_t*sin(th2);
    C1 = h_coriolis*(2*d1*d2 + d2^2);
    C2 = -h_coriolis*d1^2;
    
    G1_dyn = (m1_t*lc1_t + m2_t*l1)*g*cos(th1) + m2_t*lc2_t*g*cos(th1 + th2);
    G2_dyn = m2_t*lc2_t*g*cos(th1 + th2);
    
    tau1(i) = M11*dd1 + M12*dd2 + C1 + G1_dyn;
    tau2(i) = M12*dd1 + M22*dd2 + C2 + G2_dyn;
end

tau1_kgfcm = tau1 * (100 / 9.81);
tau2_kgfcm = tau2 * (100 / 9.81);


%% 10. Gráficas Cinemáticas y de Torque
figure(2); clf;
subplot(3,1,1); plot(t_total_vec, dtheta1, 'b', t_total_vec, dtheta2, 'r'); title('Velocidades (rad/s)'); grid on;
subplot(3,1,2); plot(t_total_vec, ddtheta1, 'b', t_total_vec, ddtheta2, 'r'); title('Aceleraciones Suavizadas (rad/s^2)'); grid on;
subplot(3,1,3); plot(t_total_vec, tau1_kgfcm, 'b', t_total_vec, tau2_kgfcm, 'r'); title('Torques (kgf·cm)'); xlabel('Tiempo [s]'); grid on;

%% 11. Animación del Robot
figure(3); clf; hold on; axis equal; grid on; axis([-0.25 0.6 -0.05 0.6]);
plot([-0.25, 0.6], [0, 0], 'k-', 'LineWidth', 3); % Mesa
plot(x_total, y_total, 'color', [0.8 0.8 0.8], 'LineWidth', 1.5); % Ruta

h_robot = plot([0, l1*cos(theta1(1)), x_real(1)], [0, l1*sin(theta1(1)), y_real(1)], '-o', 'LineWidth', 4, 'Color', '#D95319', 'MarkerFaceColor', '#0072BD');
h_rastro = plot(x_real(1), y_real(1), 'b-', 'LineWidth', 2);

x_rastro = zeros(1, num_puntos); y_rastro = zeros(1, num_puntos);
paso = 15;

disp('7. Reproduciendo animación...');
for i = 1:paso:num_puntos
    x_c = l1*cos(theta1(i)); y_c = l1*sin(theta1(i));
    x_rastro(1:i) = x_real(1:i); y_rastro(1:i) = y_real(1:i);
    
    set(h_robot, 'XData', [0, x_c, x_real(i)], 'YData', [0, y_c, y_real(i)]);
    set(h_rastro, 'XData', x_rastro(1:i), 'YData', y_rastro(1:i));
    drawnow;
end
set(h_robot, 'XData', [0, l1*cos(theta1(end)), x_real(end)], 'YData', [0, l1*sin(theta1(end)), y_real(end)]);
set(h_rastro, 'XData', x_real, 'YData', y_real); drawnow;
disp('Proceso finalizado.');