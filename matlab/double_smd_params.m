%Masses [kg]
m1 = 1.0;
m2 = 0.8;

%Spring constants [N/m]
k1 = 120;
k2 = 180;

%Damping coefficients [Ns/m]
c1 = 1.5;
c2 = 2.0;

% Continuous LTI State-space matrices
A = [0, 1, 0, 0; 
     -(k1+k2)/m1, -(c1 + c2)/m1, k2/m1, c2/m1; 
     0, 0, 0, 1; 
     k2/m2, c2/m2, -k2/m2, -c2/m2];

B = [0; 0; 0; 1/m2];

C = [1, 0, 0, 0; 
     0, 0, 1, 0];

D = [0; 0];

Ts = 1e-2; % Sampling time

save("double_smd_params.mat", "m1", "m2", ...
                              "k1", "k2", ...
                              "c1", "c2", ...
                              "A", "B", "C", "D", ...
                              "Ts");