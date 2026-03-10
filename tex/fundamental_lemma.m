clear;
clc;

n_true = 3; %prava state-space dimenzija, m = 1

% Bitno da su eigs matrice A unutar jed. kruznice - stabilnost diskretnog sust. 
A = [0.5, 0.1, 0.1; 0.1, 0.4, 0.2; -0.1, 0.2, 0.3]; % 3x3
B = [-1; 0.1; 0.7]; %3x1 = n_true x m
C = [0.1, 2, 1.5]; % 1x3
D = 0;
sys = ss(A, B, C, D, -1);

T = 50;
u = rand(1, T); %Noise input
y = lsim(sys, u, 1:T); %snimljena izlazna trajektorija

L = 22;
H_u = create_hankel(u, L); 
H_y = create_hankel(y, L);

H_data = [H_u; H_y];

% Testiranje fund. leme
x = rand(n_true, 1); % Random pocetno stanje
y_new = zeros(1, L);
u_new = rand(1, L);

% Simuliranje sustava
for k = 1:L
    y_new(k) = C*x + D*u_new(k);
    x = A*x + B*u_new(k);
end

traj_new = [u_new'; y_new']; % Ulazno - izlazna trajektorija

g = lsqminnorm(H_data, traj_new); % koristenje M-P pseudoinv.

reconstructed_traj = H_data * g;

% Prema fund. lemi traj_new mora biti u slici od H_data, tj ova razlika
% mora bit 0
disp(norm(traj_new - reconstructed_traj))