clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));

load("double_smd_params.mat");

% Constant rng seed makes the recorded output non-changing (because PRBS is
% constant) - easier for tuning the controller
% rng(1);

C = [0 0 1 0];  %Only position of mass 2
D = 0;
sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');

L_ini = 20;
L_ref = 100;

L = L_ini + L_ref;
N = 4; %Upper bound on system order, N >= n_true

T = 500; %Initial trajectory length, T >= (m + 1)(N + L) -1; m - # of inputs

[u, t] = generate_filtered_prbs(Ts, T * Ts, 80, 5, 0.03); % Filtered PRBS input

H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
shape_HNL = size(H_NL_u);

fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)); % 1 - True, 0 - False

y_c = lsim(sys_d, u, t);

sigma = 0.01;

%Make measurements different on each run of the script
rng("shuffle");

y_m = y_c + sigma * randn(size(y_c));

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y_m, L);

U_p = H_L_u(1:L_ini, :);
U_f = H_L_u(L_ini+1:L_ini + L_ref, :);
Y_p = H_L_y(1:L_ini, :);
Y_f = H_L_y(L_ini+1: L_ini + L_ref, :);

Q = 1e3;
R = 1e-3;

lambda_y = 1e4;

lambdaGValues = logspace(-8, 1, 10);
lambda_g = 5;

y_f = 1; %Step
u_f = 71.94; % Steady state input
tol = 1e-6;
u_max = 100;
y_max = 2;
deepc = RegDeePC(U_p, U_f, Y_p, Y_f, y_f, u_f, Q, R, lambda_y, lambda_g, u_max, y_max);
deepc_l = DeePC(U_p, U_f, Y_p, Y_f, y_f,u_f, Q, R, u_max, y_max);
T_sim = 250;

% Initial traj
x = repmat([0; 0; 0; 0],1, L_ini+1);
y = repmat(0, 1, L_ini);
u = repmat(0, 1, L_ini);

for k=1:T_sim
    [u_i, info] = deepc.solve(u(:, end-L_ini+1:end), y(:, end-L_ini+1:end));
    if info.errorcode ~= 0
        error("Problem occurred: %s\n", yalmiperror(info.errorcode));
    end

    u = [u, u_i];

    y_clean = C * x(:, end);

    y_noisy = y_clean + sigma * randn;

    y = [y, y_noisy];    
    x = [x, sys_d.A * x(:, end) + sys_d.B * u(:, end)];

    % if norm(y(:, end) - y_f, inf) < tol
    %     break;
    % end
end

x = repmat([0; 0; 0; 0],1, L_ini+1);
y_l = repmat(0, 1, L_ini);
u_l = repmat(0, 1, L_ini);

for k=1:T_sim
    [u_i, info] = deepc_l.solve(u_l(:, end-L_ini+1:end), y_l(:, end-L_ini+1:end));
    if info.errorcode ~= 0
        error("Problem occurred: %s\n", yalmiperror(info.errorcode));
    end

    u_l = [u_l, u_i];

    y_clean = C * x(:, end);

    y_noisy = y_clean + sigma * randn;

    y_l = [y_l, y_noisy];    
    x = [x, sys_d.A * x(:, end) + sys_d.B * u_l(:, end)];

    % if norm(y(:, end) - y_f, inf) < tol
    %     break;
    % end
end

subplot(2, 1, 1);
plot((1:length(y))*Ts, y, 'LineWidth',2.0, 'Color',[0 0 1], 'LineStyle','-');
hold on;
plot((1:length(y))*Ts, y_l, 'LineWidth',2.0, 'Color',[0.5 0 0], 'LineStyle','-');
hold off;
xlim([L_ini*Ts T_sim*Ts]);
ylim([0 2]);
yline(1, 'LineWidth',1.0, 'Color',[0 0.5 0], 'LineStyle','--');
xlabel("Time [s]");
ylabel("Mass 2 position [m]");
legend("Regularized DeePC", "DeePC");
grid on;
ax = gca;
ax.FontSize = 14;

subplot(2, 1, 2);
plot((1:length(y))*Ts, u, 'LineWidth',2.0, 'Color',[1 0 1], 'LineStyle','-');
hold on;
plot((1:length(y))*Ts, u_l, 'LineWidth',2.0, 'Color',[0 0.5 0], 'LineStyle','-');
hold off;
xlim([L_ini*Ts T_sim*Ts]);
ylim([-120 120]);
yline(0, 'LineWidth',1.0, 'Color',[0.5 0 0], 'LineStyle','--');
xlabel("Time [s]");
ylabel("Input force [N]");
legend("Regularized DeePC", "DeePC");
grid on;
ax = gca;
ax.FontSize = 14;

disp(y(end));