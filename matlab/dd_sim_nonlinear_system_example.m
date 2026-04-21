%Comparison between data-driven and model based simulations of driven pendulum system
clear; clc;

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

[u, t] = generate_filtered_prbs(Ts, T * Ts, 1, 3, 0.01);

for k = 1:T-1
    y(k+1, :) = pendulum_rk4(y(k, :), u(k), Ts, m, L, b, g); 
end

plot(t, y(:, 1), '');
ylim([-pi pi]);
xlim([0 50]);
grid on;

function ynext = pendulum_rk4(y, u, Ts, m, L, b, g)
    f = @(x,u) [ 
        x(2),...
       -(b/(m*L^2))*x(2) - (g/L)*sin(x(1)) + u/(m*L^2)
    ];

    k1 = f(y, u);
    k2 = f(y + 0.5*Ts*k1, u);
    k3 = f(y + 0.5*Ts*k2, u);
    k4 = f(y + Ts*k3, u);

    ynext = y + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end
