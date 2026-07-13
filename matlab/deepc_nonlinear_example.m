clear; clc; addpath(fullfile(fileparts(mfilename("fullpath")), "custom_functions"));
load("double_smd_params.mat");

C = [0 0 1 0];
k23 = 24;

Ts = 0.01;
u = 100;
x = zeros(4, 1);
x_l = zeros(4, 1);

y = zeros(100, 1);
y_l = zeros(100, 1);

for i = 1:100
    y(i) = C*x;
    y_l(i) = C*x_l;
    x = nonlinear_double_smd(x, u, Ts, A, B, m1, m2, k23);
    x_l = nonlinear_double_smd(x_l, u, Ts, A, B, m1, m2, 0);
end

plot((1:100)*Ts, y_l);
hold on;
plot((1:100)*Ts, y);
hold off;
legend("Linear", "Nonlinear");