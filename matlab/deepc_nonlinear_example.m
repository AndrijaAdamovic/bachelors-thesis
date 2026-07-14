clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));
load("double_smd_params.mat");

% Constant rng seed makes the recorded output non-changing (because PRBS is
% constant) - easier for tuning the controller
rng(2);

C = [0 0 1 0];  %Only position of mass 2
k23 = 40;

L_ini = 20;
L_ref = 100;

L = L_ini + L_ref;
N = 4; %Upper bound on system order, N >= n_true

T = 500; %Initial trajectory length, T >= (m + 1)(N + L) -1; m - # of inputs

[u, t] = generate_filtered_prbs(Ts, T * Ts, 100, 5, 0.03);

H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
shape_HNL = size(H_NL_u);

fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)); % 1 - True, 0 - False

sigma = 0.01;

x = zeros(4, 1);
y = zeros(500, 1);

for i = 1:T
    y(i) = C*x;
    x = nonlinear_double_smd(x, u(i), Ts, A, B, m1, m2, k23);
end

%Make measurements different on each run of the script
rng("shuffle");

y = y + sigma .* randn(size(y));

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

U_p = H_L_u(1:L_ini, :);
U_f = H_L_u(L_ini+1:L_ini + L_ref, :);
Y_p = H_L_y(1:L_ini, :);
Y_f = H_L_y(L_ini+1: L_ini + L_ref, :);

Q = 1e3;
R = 1e-3;

lambda_y = 1e7;

lambdaGValues = logspace(-8, 1, 10);
lambda_g = 50;

y_f = 1; %Step
u_f = 72; % Steady state input
tol = 1e-6;
u_max = 100;
y_max = 2;
deepc = RegDeePC(U_p, U_f, Y_p, Y_f, y_f, u_f, Q, R, lambda_y, lambda_g, u_max, y_max);

T_sim = 200;

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
    x = [x, nonlinear_double_smd(x(:, end), u(:, end), Ts, A, B, m1, m2, k23)];

    % if norm(y(:, end) - y_f, inf) < tol
    %     break;
    % end
end

subplot(2, 1, 1);
plot((1:length(y))*Ts, y, 'LineWidth',2.0, 'Color',[0 0 1], 'LineStyle','-');
xlim([L_ini*Ts T_sim*Ts]);
ylim([0 2]);
yline(1, 'LineWidth',1.0, 'Color',[0 0.5 0], 'LineStyle','--');
ylabel("Mass 2 position [m]");
grid on;
ax = gca;
ax.FontSize = 14;

subplot(2, 1, 2);
plot((1:length(y))*Ts, u, 'LineWidth',2.0, 'Color',[1 0 1], 'LineStyle','-');
xlim([L_ini*Ts T_sim*Ts]);
ylim([-120 120]);
yline(0, 'LineWidth',1.0, 'Color',[0.5 0 0], 'LineStyle','--');
ylabel("Input force [N]");
xlabel("Time [s]");
grid on;
ax = gca;
ax.FontSize = 14;