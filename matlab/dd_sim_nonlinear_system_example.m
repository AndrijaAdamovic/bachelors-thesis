%Comparison between data-driven and model based simulations of driven pendulum system
clear; clc; addpath('custom_functions');

%Pendulum parameters
m = 0.2;
L = 0.5;
b = 0.05;
g = 9.81;

%Sampling time
Ts = 1e-2;


T = 50 / Ts;

%State
y = zeros(T, 2);
y(1, :) = [pi/2; 0]; % Inital state

[u, t] = generate_filtered_prbs(Ts, T * Ts, 1, 2, 0.01);

for k = 1:T-1
    y(k+1, :) = pendulum_rk4(y(k, :), u(k), Ts, m, L, b, g); 
end

plot(t, y(:, 1), '');
ylim([-pi pi]);
xlim([0 50]);
grid on;