% Data-driven simulation of system from "fl_linear_system_example_1.m"
clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));

load("double_smd_params.mat");

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');

L_ini = 5 / Ts;
L_ref = 10 / Ts;

L = L_ini + L_ref;
N = 4; %Upper bound on system order, N >= n_true

T = 50 / Ts; %Initial trajectory length, T >= (m + 1)(N + L) -1; m - # of inputs

[u, t] = generate_filtered_prbs(Ts, T * Ts, 2, 5, 0.03); % Filtered PRBS input

H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
shape_HNL = size(H_NL_u);

fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)); % 1 - True, 0 - False

y = lsim(sys_d, u, t);

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

U_p = H_L_u(1:L_ini, :);
U_f = H_L_u(L_ini+1:L_ini + L_ref, :);
Y_p = H_L_y(1:2*L_ini, :);
Y_f = H_L_y(2*L_ini+1: 2*L_ini + 2*L_ref, :);

M = [U_p; Y_p; U_f];

t_s = (0:L-1) * Ts;
u_s = 10 * (sin(2 * pi * 1 * t_s) + (1/3)*sin(2 * pi * 3 * t_s));
y_s = lsim(sys_d, u_s, t_s);
y_s_vec = reshape(y_s.', [], 1);

b = [u_s(1:L_ini)'; y_s_vec(1:2*L_ini); u_s(L_ini+1:L)'];

g = lsqminnorm(M, b);

y_f_vec = Y_f * g;
y_f = reshape(y_f_vec, 2, []).';

plot(t_s, y_s(:, 1));
hold on;
plot(t_s(L_ini+1:end), y_f(:, 1));
hold off;

legend("Model sim", "Data-based")  