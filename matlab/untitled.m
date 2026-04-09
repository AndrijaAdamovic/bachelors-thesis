clear; clc;

%Pendulum parameters
m = 0.2;
L = 0.5;
b = 0.05;
g = 9.81;

Ts = 1e-2;
N = 10 / Ts;
t = (0:N)*Ts;
u = t * 0;

x = zeros(2, N+1);
x(:, 1) = [pi/2; 0];

for k = 1:N
    x(:, k+1) = pendulum_rk4(x(:,k), u(k), Ts, m, L, b, g); 
end

plot(t, x(1, :));

function xnext = pendulum_rk4(x, u, Ts, m, L, b, g)
    f = @(x,u) [ 
        x(2);
       -(b/(m*L^2))*x(2) - (g/L)*sin(x(1)) + u/(m*L^2)
    ];

    k1 = f(x, u);
    k2 = f(x + 0.5*Ts*k1, u);
    k3 = f(x + 0.5*Ts*k2, u);
    k4 = f(x + Ts*k3, u);

    xnext = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end
