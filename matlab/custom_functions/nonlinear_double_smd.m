function xnext = nonlinear_double_smd(x, u, Ts, A, B, m1, m2, k23)
    f = @(x, u) A*x + B*u + [0; k23/m1; 0; k23/m2]*(x(3) - x(1))^3;

    k1 = f(x, u);
    k2 = f(x + 0.5*Ts*k1, u);
    k3 = f(x + 0.5*Ts*k2, u);
    k4 = f(x + Ts*k3, u);

    xnext = x + (Ts/6)*(k1 + 2*k2 + 2*k3 + k4);
end

