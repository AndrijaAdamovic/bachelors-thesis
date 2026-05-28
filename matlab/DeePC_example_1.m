% DeePC example for the double spring mass damper system
clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));

load("double_smd_params.mat");

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');

L_ini = round(0.2 / Ts);
L_ref = round(0.4 / Ts);
L = L_ini + L_ref;

T = 5 / Ts; %Initial trajectory length, T >= (m + 1)(N + L) -1; m - # of inputs

N = 4; %Upper bound on system order, N >= n_true

[u, t] = generate_filtered_prbs(Ts, T * Ts, 2, 5, 0.03); % Filtered PRBS input

% Test for persistency of excitation
% H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
% shape_HNL = size(H_NL_u);
% 
% fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)); % 1 - True, 0 - False

y = lsim(sys_d, u, t); 

m = size(u, 2); % number of inputs
p = size(y, 2); % number of outputs

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

U_p = H_L_u(1 : m*L_ini, :);
U_f = H_L_u(m*L_ini + 1 : m*L, :);
Y_p = H_L_y(1 : p*L_ini, :);
Y_f = H_L_y(p*L_ini + 1: p * L, :);

% DeePC Sim
yalmip('clear');

n_g = size(U_p, 2); % number of columns of hankel matrix i.e. length of g

g = sdpvar(n_g, 1);

u_pred = U_f * g;
y_pred = Y_f * g;

Qy = eye(p);
R = 1e-3 * eye(m);

Qbar = kron(eye(L_ref), Qy);
Rbar = kron(eye(L_ref), Ru);

u_max = 2; % Input constraint



