%Comparison between data-driven and model based simulations of driven pendulum system
clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));

%Pendulum parameters
m = 0.2;
l = 0.5;
b = 0.05;
g = 9.81;

%Sampling time
Ts = 1e-2;

T = 50 / Ts;

%State - both angle and angular velocity 
y = zeros(T, 2);
y(1, :) = [0; 0]; % Inital state

[u, t] = generate_filtered_prbs(Ts, T * Ts, 0.1, 2, 0.05);

for k = 1:T-1
    y(k+1, :) = pendulum_rk4(y(k, :), u(k), Ts, m, l, b, g); 
end

L_ini = 2 / Ts;
L_ref = 20 / Ts;

L = L_ini + L_ref;
N = 4; % Upper bound of system order,linearized system has order 2  

fprintf("Sufficiently long trajectory, T >= (m + 1)(N + L) -1? : %d \n", T >= (m + 1)*(N + L) - 1)
 
H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

U_p = H_L_u(1:L_ini, :);
U_f = H_L_u(L_ini+1:L_ini + L_ref, :);
Y_p = H_L_y(1:2*L_ini, :);
Y_f = H_L_y(2*L_ini+1: 2*L_ini + 2*L_ref, :);

M = [U_p; Y_p; U_f];

%Simulate system from initial state of [pi/2 0], and no input
t_s = (0:L-1) * Ts;
u_s = zeros(1, L);
y_s = zeros(L, 2);
y_s(1, :) = [pi/10 0];

for k = 1:L-1
    y_s(k+1, :) = pendulum_rk4(y_s(k, :), 0, Ts, m, l, b, g); 
end

y_s_vec = reshape(y_s.', [], 1); % Vectorize data

b_d = [u_s(1:L_ini)'; y_s_vec(1:2*L_ini); u_s(L_ini+1:L)'];

g_vec = lsqminnorm(M, b_d);

y_f_vec = Y_f * g_vec;
y_f = reshape(y_f_vec, 2, []).';

plot(t_s, y_s(:, 1), '');
ylim([-pi pi]);
xlim([0 t_s(end)]);
grid on;
hold on;
plot(t_s(L_ini+1:end), y_f(:, 1));

% plot(t, y(:, 1), '');
% ylim([-pi pi]);
% xlim([0 t(end)]);
% yticks(-pi:pi/20:pi)
% grid on;