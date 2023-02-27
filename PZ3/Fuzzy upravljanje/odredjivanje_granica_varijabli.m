%clc 
clear
close all
%% Parametri modela
tfin=60;     % sec
V = 4; %l
SF = 10; %g/l
Y = 0.5;
mi_max = 1;
K1 = 0.03; %g/l
K2 = 0.5; %l/g
%% Nominalni režim, ranije izračunato
Se = 0.2187;
Xe = 4.8907;
Fe = 3.2089;

%% Limit upravljanja 
F_upper_limit = Fe+1e-3*Fe;
F_lower_limit = 0;
%%  Početno stanje varijabli: F_lower_limit
x10 = Y*SF;
x20 = 0; 
%%
sim('granice_nelin_sim')
%% Odzivi u eksperimentu sa stepom upravljanja Fmin->Fmax
figure;
plot(t_out, S_out); xlabel('t [h]'); ylabel('S [g/l]');
figure;
plot(t_out, derivative_S_out); xlabel('t [h]'); ylabel('dS/dt [g/l/h]');

figure;
plot(t_out, derivative_F_out); xlabel('t [h]'); ylabel('dF/dt [l/h/h]');

%%
S_upper_limit = max(S_out)
S_lower_limit = min(S_out)

derivative_S_upper_limit = max(derivative_S_out)
derivative_S_lower_limit = min(derivative_S_out)
%%

e_upper_limit = Se - S_lower_limit
e_lower_limit = Se - S_upper_limit

derivative_e_upper_limit = -derivative_S_lower_limit;
derivative_e_lower_limit = -derivative_S_upper_limit;


%%
derivative_F_upper_limit = 3.2121e+04;
derivative_F_lower_limit = -3.2121e+04;