% -------------------------------------------------------------
% Sistema discreto descrito por diagrama de bloques
% -------------------------------------------------------------

clear; clc; close all;
Ts = 1;  % tiempo de muestreo (unidad genérica)

%% a) Modelo original en variables de estado
A = [0  -0.25;
     1   0    ];
B = [1; 0];
C = [3  -4;
     2  -2;
     0  -1];
D = zeros(3,1);
sys_orig = ss(A, B, C, D, Ts);
disp('a) Matrices del modelo original:');
disp('A = '); disp(A);
disp('B = '); disp(B);
disp('C = '); disp(C);
disp('D = '); disp(D);

%% b) Matriz de transferencia H(z)
z = tf('z', Ts);
H1 = (3*z - 4) / (z^2 + 0.25);
H2 = (2*z - 2) / (z^2 + 0.25);
H3 = -1 / (z^2 + 0.25);
H = [H1; H2; H3];
disp('b) Matriz de transferencia H(z):');
disp(H);

%% c) Nuevo modelo a partir de H(z)
sys_new = ss(H);
A2 = sys_new.A;
B2 = sys_new.B;
C2 = sys_new.C;
D2 = sys_new.D;
disp('c) Matrices del nuevo modelo:');
disp('A2 = '); disp(A2);
disp('B2 = '); disp(B2);
disp('C2 = '); disp(C2);
disp('D2 = '); disp(D2);

%% d) Verificación de equivalencia (respuesta al escalón)
figure;
step(sys_orig, 'b--', sys_new, 'r-');
legend('Original', 'Nuevo', 'Location', 'best');
title('d) Comparación de respuestas al escalón');
grid on;

% Transformación de semejanza para mapear condiciones iniciales
Co_orig = ctrb(A, B);
Co_new  = ctrb(A2, B2);
T = Co_new / Co_orig;   % x_nuevo = T * x_orig

%% e) Simulaciones con el modelo nuevo
N = 15;                     % número de muestras
t = (0:N-1)';               % vector de tiempo

% Condiciones iniciales originales y su equivalente en el nuevo modelo
x0_orig = [1; -1];
x0_new = T * x0_orig;

% e.1) Respuesta con entrada cero y CI dadas
[y0, t0, x0] = initial(sys_new, x0_new, N-1);
figure;
subplot(3,2,1); stairs(t0, y0(:,1)); ylabel('y_1'); title('e.1) Salida y_1 (entrada cero)'); grid on;
subplot(3,2,3); stairs(t0, y0(:,2)); ylabel('y_2'); title('Salida y_2'); grid on;
subplot(3,2,5); stairs(t0, y0(:,3)); ylabel('y_3'); xlabel('n'); grid on;
subplot(3,2,2); stairs(t0, x0(:,1)); ylabel('x_{new,1}'); title('Estado x_{new,1}'); grid on;
subplot(3,2,4); stairs(t0, x0(:,2)); ylabel('x_{new,2}'); xlabel('n'); grid on;

% e.2) Respuesta al escalón sin CI
figure;
step(sys_new, N-1);
title('e.2) Respuesta al escalón (condiciones iniciales nulas)');
grid on;

% e.3) Respuesta al escalón con CI dadas
u_step = ones(N, 1);
[y_step, t_step, x_step] = lsim(sys_new, u_step, t, x0_new);
figure;
subplot(3,2,1); stairs(t_step, y_step(:,1)); ylabel('y_1'); title('e.3) Salidas con escalón + CI'); grid on;
subplot(3,2,3); stairs(t_step, y_step(:,2)); ylabel('y_2'); grid on;
subplot(3,2,5); stairs(t_step, y_step(:,3)); ylabel('y_3'); xlabel('n'); grid on;
subplot(3,2,2); stairs(t_step, x_step(:,1)); ylabel('x_{new,1}'); title('Estados'); grid on;
subplot(3,2,4); stairs(t_step, x_step(:,2)); ylabel('x_{new,2}'); xlabel('n'); grid on;