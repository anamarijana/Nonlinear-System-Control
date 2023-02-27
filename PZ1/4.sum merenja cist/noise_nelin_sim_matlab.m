clc;
clear;
close all;
%% Nominalni parametri modela

tfin=15;     % sec

V = 4; %l
SF = 10; %g/l
Y = 0.5;
mi_max = 1;
K1 = 0.03; %g/l
K2 = 0.5; %l/g
%% Nominalni režim (F*X)max
syms s_Se

mi =        mi_max*s_Se/(K2*s_Se^2 + s_Se + K1);
fx_max = s_Se*(SF-s_Se)/(K2*s_Se^2 + s_Se + K1); %funkcija koju maksimizujem
dif_fx_max = diff(fx_max,s_Se); 

Ses = double(solve(dif_fx_max, s_Se));
Se = Ses(2);

s_Se = Se;
mi_Se = eval(mi);

Fe = double(mi_Se*V);
Xe = Y*(SF-Se);

%% Ponašanje oko nominalnog režima

S_initial = Se+0.1*Se; %proizvoljna tačka dovoljno blizu
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);

x10 = X_initial;
x20 = S_initial;


F = Fe;

%% Uvođenje šuma merenja

std_noise1 = 1e-3*Se;
std_noise2 = 1e-4*Se;
std_noise3 = 1e-5*Se;

%% Pokretanje simulacije
sim('noise_nelin_sim.slx');

%%
%JEDINO NA ŠTA UTIČE ŠUM MERENJA U OTVORENOJ SPREZI
figure;
hold all;
plot(t_out,S_out,'r')%,'LineWidth',2)
plot(t_out,S_out1,'g')%,'LineWidth',2)
plot(t_out,S_out2,'b')%,'LineWidth',2)
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
legend('\sigma_{noise} = 1e-3Se','\sigma_{noise} = 1e-4Se','\sigma_{noise} = 1e-5Se')
title('Koncentracija supstrata na izlazu')