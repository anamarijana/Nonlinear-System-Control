clc;
clear;
close all;
%% Nominalni parametri modela

tfin=30;     % sec

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

% x10 = X_initial;
% x20 = S_initial;
%% Limit upravljanja 
F_upper_limit = 1e-3*Fe;
F_lower_limit = 0;

%%

x10 = Xe;
x20 = Se;
F = Fe;

F_step1 = 0;
F_step2 = Fe;
tF_step1 = 2;
tF_step2 = 5;




%% Pokretanje simulacije
sim('nelin_sim.slx');


%% Procena propusnog opsega
response = -S_out+S_out(1);
yend = response(end);
idx = find(response<=0.63*abs(yend));
max_idx = idx(end);
time_constant = t_out(max_idx);
disp('Vremenska konstanta sistema:')
disp(time_constant)
w0 = 1./time_constant;
disp('Propusni opseg:')
disp(w0)

%% procenjeno sa grafika
S_upper_limit = 0.2325;
S_lowe_limit = 0;
%% Rezultati 

%PROMENLJIVE STANJA U VREMENU

figure;
subplot(2,1,1)
plot(t_out,x1_out,'k')%,'LineWidth',2)
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
grid

subplot(2,1,2)
plot(t_out,x2_out,'k')%,'LineWidth',2)
xlabel('vreme [h]')
ylabel('x2 = S [g/l]')
sgtitle('Promenljive stanja')
grid

%% 

%FAZNI DIJAGRAM

figure;
plot(x1_out,x2_out,'k')%,'LineWidth',2)
xlabel('x1 = X [g/l]')
ylabel('x2 = S [g/l]')
title('Fazni digram')
grid
%%

%UPRAVLJANJE I IZLAZ

figure;
subplot(2,1,1)
plot(t_out,F_out,'k')%,'LineWidth',2)
grid
title('Zapremisnki protok kroz reaktor')
ylabel('u = F [l/h]')
xlabel('vreme [h]')

subplot(2,1,2)
plot(t_out,S_out,'k')%,'LineWidth',2)
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
title('Koncentracija supstrata na izlazu')
% %
% figure;
% plot(t_out,response,'k')%,'LineWidth',2)
% grid
% ylabel('y = S [g/l]')
% xlabel('vreme [h]')
% title('Koncentracija supstrata na izlazu')
% 
% 
% figure;
% subplot(2,1,1)
% plot(t_out,F_lin_out,'k','LineWidth',2)
% grid
% title('Linearni sistem: Protok supstrata u reaktor')
% ylabel('F [l/h]')
% xlabel('vreme [h]')
% 
% subplot(2,1,2)
% plot(t_out,S_lin_out,'k','LineWidth',2)
% grid
% ylabel('S [g/l]')
% xlabel('vreme [h]')
% title('Linearni sistem: Koncentracija supstrata na izlazu')

