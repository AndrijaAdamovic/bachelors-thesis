% Dual mass spring damper system example for demonstrating the Fundamental
% lemma
clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));

load("double_smd_params.mat");

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');

N = 4; %Upper bound on system order, N >= n_true
L = round(12 / Ts); %Length of reconstructed trajectory, lasting 12 seconds

T = round(30 / Ts); %Initial trajectory length, T >= (m + 1)(N + L) -1; m - # of inputs

[u, t] = generate_filtered_prbs(Ts, T * Ts, 2, 5, 0.03); % Filtered PRBS input

%Test persistency of excitation
H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
shape_HNL = size(H_NL_u);

fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)); % 1 - True, 0 - False

%Input PRBS
y = lsim(sys_d, u, t);

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

H_data = [H_L_u; H_L_y];

%Step input
u_step = ones(L, 1);
y_step = lsim(sys_d, u_step);

y_step_vec = reshape(y_step.', [], 1); % Vectorization

step_trajectory = [u_step; y_step_vec];

%Solve for vector g using least squares
g_step = lsqminnorm(H_data, step_trajectory);

recon_step_trajectory = H_data * g_step;

fprintf("||step_traj - reconstructed_traj|| = %d \n", norm(step_trajectory - recon_step_trajectory))
recon_u_step = recon_step_trajectory(1:length(u_step));
recon_y_step = reshape(recon_step_trajectory(length(u_step)+1:end), 2, []).';


% plot((1:L), y_step(:, 1), 'LineWidth', 2.0);
% hold on;
% plot((1:L), recon_y_step(:, 1), 'LineStyle', '--', 'LineWidth', 2.0);
% hold off;

% PRBS reconstruction
u_prbs = generate_filtered_prbs(Ts, L * Ts, 4, 10, 0.02); 
y_prbs = lsim(sys_d, u_prbs);

y_prbs_vec = reshape(y_prbs.', [], 1);

prbs_trajectory = [u_prbs; y_prbs_vec];

g_prbs = lsqminnorm(H_data, prbs_trajectory);
recon_prbs_trajectory = H_data * g_prbs;
fprintf("||prbs_traj - reconstructed_traj|| = %d \n", norm(prbs_trajectory - recon_prbs_trajectory))

recon_u_prbs = recon_prbs_trajectory(1:length(u_prbs));
recon_y_prbs = reshape(recon_prbs_trajectory(length(u_prbs)+1:end), 2, []).';

% plot((1:L), y_prbs(:, 2), 'LineWidth', 2.0);
% hold on;
% plot((1:L), recon_y_prbs(:, 2), 'LineStyle', '--', 'LineWidth', 2.0);
% hold off;