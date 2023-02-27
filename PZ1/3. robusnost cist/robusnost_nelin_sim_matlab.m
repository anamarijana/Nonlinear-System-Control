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

%% Ispitivanje robusnosti
Y1 = Y - 0.2*Y;
Y2 = Y - 0.1*Y;
Y3 = Y + 0.1*Y;
Y4 = Y + 0.2*Y;
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

S_initial = Se + 0.1*Se; %proizvoljna tačka dovoljno blizu
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);

x10 = X_initial;
x20 = S_initial;

F = Fe;
%% Pokretanje simulacije
open('robusnost_nelin_sim.slx')
sim('robusnost_nelin_sim.slx');
%% Rezultati
%PROMENLJIVE STANJA U VREMENU
figure;
subplot(2,1,1)
hold all;
plot(t_out,x1_out)
plot(t_out,x1_out1)
plot(t_out,x1_out2)
plot(t_out,x1_out3)
plot(t_out,x1_out4)
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
grid

subplot(2,1,2)
hold all;
plot(t_out,x2_out)
plot(t_out,x2_out1)
plot(t_out,x2_out2)
plot(t_out,x2_out3)
plot(t_out,x2_out4)
xlabel('vreme [h]')
ylabel('x2 = S [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
sgtitle('Promenljive stanja')
grid

%% 

%FAZNI DIJAGRAM

figure;
hold all;
plot(x1_out,x2_out)
plot(x1_out1,x2_out1)
plot(x1_out2,x2_out2)
plot(x1_out3,x2_out3)
plot(x1_out4,x2_out4)
xlabel('x1 = X [g/l]')
ylabel('x2 = S [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
title('Fazni digram')
grid
%%

%UPRAVLJANJE I IZLAZ

figure;
subplot(2,1,1)
hold all;
plot(t_out,F_out)
plot(t_out,F_out1)
plot(t_out,F_out2)
plot(t_out,F_out3)
plot(t_out,F_out4)
grid
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
title('Zapremisnki protok kroz reaktor')
ylabel('u = F [l/h]')
xlabel('vreme [h]')

subplot(2,1,2)
hold all;
plot(t_out,S_out)
plot(t_out,S_out1)
plot(t_out,S_out2)
plot(t_out,S_out3)
plot(t_out,S_out4)
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
title('Koncentracija supstrata na izlazu')


%% Procena propusnog opsega 0.8 0.9 1 
% 
%nominalni
response1 = -S_out+S_out(1);
yend = response1(end);
idx = find(response1<=0.63*abs(yend));
max_idx = idx(end);
time_constant = t_out(max_idx);
w0_estimated = 1./time_constant;
disp('Propusni opseg regulacije, nominalni Y:')
disp(w0_estimated)
%%
%0.8
response1 = S_out1;
yend = response1(end);
idx = find(response1<=0.63*abs(yend));
max_idx1 = idx(end);
time_constant1 = t_out(max_idx1);
w0_estimated1 = 1./time_constant1;
disp('Propusni opseg regulacije, 0.8 Y nominalnog:')
disp(w0_estimated1)
%%
%0.9
response2 = S_out2;
yend2 = response2(end);
idx = find(response2<=0.63*abs(yend2));
max_idx2 = idx(end);
time_constant2 = t_out(max_idx2);
w0_estimated2 = 1./time_constant2;
disp('Propusni opseg regulacije,0.9 Y nominalnog:')
disp(w0_estimated2)

%% samo koji valjaju 

figure;
hold all;
plot(t_out,S_out)
plot(t_out,S_out1)
plot(t_out,S_out2)
%plot(t_out(max_idx),S_out(max_idx), 'r.')
%plot(t_out(max_idx2),S_out2(max_idx2), 'r.')
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
title('Koncentracija supstrata na izlazu')

figure;
hold all;
plot(t_out,x1_out)
plot(t_out,x1_out1)
plot(t_out,x1_out2)
ylim([3 5])
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
title('Koncentracija biomase na izlazu')
grid
