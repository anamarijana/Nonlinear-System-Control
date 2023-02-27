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
%% Poremećaj u toku nominalnog režima

x10 = Xe;
x20 = Se;
F = Fe;

%% Uvođenje poremećaja
disturb_X1 = 0.2*Xe;
disturb_X2 = -0.2*Xe;
t_disturb = 2;
%% Pokretanje simulacije
open('poremecaj_nelin_sim.slx')
sim('poremecaj_nelin_sim.slx');
%% Rezultati
%PROMENLJIVE STANJA U VREMENU
figure;
hold all;
plot(t_out,x1_out)
plot(t_out,x1_out1)
plot(t_out,x1_out2)
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
legend('d = 0','d = 0.2Xe','d = -0.2Xe')
title('Koncentracija biomase na izlazu')
grid

% subplot(2,1,2)
% hold all;
% plot(t_out,x2_out)
% plot(t_out,x2_out1)
% plot(t_out,x2_out2)
% xlabel('vreme [h]')
% ylabel('x2 = S [g/l]')
% legend('d = 0','d = 0.2Xe','d = -0.2Xe')
% sgtitle('Promenljive stanja')
% grid

%%

%UPRAVLJANJE I IZLAZ

% figure;
% subplot(2,1,1)
% hold all;
% plot(t_out,F_out)
% plot(t_out,F_out1)
% plot(t_out,F_out2)
% grid
% legend('d = 0','d = 0.2Xe','d = -0.2Xe')
% title('Zapremisnki protok kroz reaktor')
% ylabel('u = F [l/h]')
% xlabel('vreme [h]')

figure;
hold all;
plot(t_out,S_out)
plot(t_out,S_out1)
plot(t_out,S_out2)
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
legend('d = 0','d = 0.2Xe','d = -0.2Xe')
title('Koncentracija supstrata na izlazu')


%% d = 0.2Xe

response1 = S_out1;
yend = response1(end);
idx = find(response1<=0.63*abs(yend));
max_idx1 = idx(end);
time_constant1 = t_out(max_idx1)-t_disturb;
w0_estimated1 = 1./time_constant1;
disp('Propusni opseg')
disp(w0_estimated1)
% figure;
% hold all;
% plot(t_out,response1)
% plot(t_out(max_idx1),response1(max_idx1),'r.')