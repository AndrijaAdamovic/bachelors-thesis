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
