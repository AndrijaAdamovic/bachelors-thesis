clear; clc; addpath('custom_functions');

A = [-2 -6.25; 4 0];
B = [2; 0];
C = [0 3.125];
D = 0;
Ts = 1e-2;

sys_c = ss(A, B, C, D);
sys_d = c2d(sys_c, Ts, 'zoh');

T = round(20 / Ts);
N = 2;
[u, t] = generate_filtered_prbs(Ts, T * Ts, 2, 5, 0.01);

L_ini = 0.5 / Ts;
L_ref = 3 / Ts;
L = L_ini + L_ref;

H_NL_u = create_hankel(u, N + L); %Hankel matrix of depth N + L
shape_HNL = size(H_NL_u);

fprintf("Persistently exciting of order N + L: %d \n", shape_HNL(1) == rank(H_NL_u)) % 1 - True, 0 - False

y = lsim(sys_d, u, t);

H_L_u = create_hankel(u, L);
H_L_y = create_hankel(y, L);

U_p = H_L_u(1:L_ini, :);
U_f = H_L_u(L_ini+1:L_ini + L_ref, :);
Y_p = H_L_y(1:L_ini, :);
Y_f = H_L_y(L_ini+1:L_ini + L_ref, :);

t_s = (0:Ts:(L-1)*Ts);
u_s = ones(size(t_s));
% u_s = 10*sin(2 * pi * 3 * t_s);
y_s = lsim(sys_d, u_s, t_s);

b = [u_s(1:L_ini)'; y_s(1:L_ini); u_s(L_ini+1:L)'];
Q = [U_p; Y_p; U_f];

g = lsqminnorm(Q, b);

y_f = Y_f * g;

plot(t_s, y_s, 'Color','r', 'LineWidth',1);
hold on;
plot(t_s(L_ini+1:L), y_f, 'Color','b','LineWidth',1);
hold off;