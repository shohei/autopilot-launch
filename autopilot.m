clc; close all; clear all;

lyy = 2.186e8;
m = 38901;
Tc = 2.361e6;
V = 1347;
Cn_alpha = 0.1465;
g = 26.10;
N_alpha = 686819;
M_alpha = 0.3807;
M_delta = 0.5726;
x_cg = 53.19;
x_cp = 121.2;
F = Tc;

Mach = 1.4;
h = 34000;
S = 116.2;
Fbase = 1000;
Ca = 2.4;
D = Ca*680*S - Fbase;
Drag = 7.15*D;

A_m = [0 1 0; 
    M_alpha 0 M_alpha/V; 
    -(F-Drag+N_alpha)/m 0 -N_alpha/(m*V)];
B_m = [0; M_delta; Tc/m];
C_m = diag([1 1 1]);
D_m = [0;0;0];
pitch_ss = ss(A_m,B_m,C_m,D_m);


cvector = {'bo','ro','go'}';
R_vector = [0.1 5 10];

figure; hold on;

for k=1:1
   R_matrix_drift = R_vector(k);
   Q_matrix_drift = [1 0 1/V;0 0 0; 1/V 0 1/V^2];   
   [K S e] = lqr(pitch_ss, Q_matrix_drift, R_matrix_drift);
   
   for i=1:100000
      e_val(:,i) = eig(A_m-B_m*K*i/10000);       
   end
   
   plot(real(e_val(1,:)),imag(e_val(1,:)),cvector{k});
   plot(real(e_val(2,:)),imag(e_val(2,:)),cvector{k});
   plot(real(e_val(3,:)),imag(e_val(3,:)),cvector{k});   
   grid;
    
end

xlim([-2 1]);
legend('R=0.1');

%LQR gains
K_1 = K(1);
K_2 = K(2);
K_3 = K(3);







































