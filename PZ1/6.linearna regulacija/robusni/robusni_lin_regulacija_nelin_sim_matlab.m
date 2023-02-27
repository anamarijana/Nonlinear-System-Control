clc;
clear;
close all;
%% Parametri modela

tfin=60;     % sec
V = 4; %l
SF = 10; %g/l
Ynom = 0.5;
Y = Ynom;
mi_max = 1;
K1 = 0.03; %g/l
K2 = 0.5; %l/g

%% Ispitivanje robusnosti
Y1 = Y - 0.2*Y;
Y2 = Y - 0.1*Y;
Y3 = Y + 0.1*Y;
Y4 = Y + 0.2*Y;
%% Nominalni režim, ranije izračunato

Se = 0.218662692463450;
Xe = 4.89066865376828;
Fe = 3.20891060125344;
%%
s = tf('s');
G = 2.44533432688414/(s+0.802227650313360);

%% Limit upravljanja 
F_upper_limit = 1e-3*Fe;
F_lower_limit = 0;

% procenjeno sa grafika
S_upper_limit = 0.2325;
S_lowe_limit = 0;

%% tranzijent u nominalni

S_initial = S_upper_limit; %proizvoljna tačka dovoljno blizu
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);
% 
% x10 = X_initial;
% x20 = S_initial;


%% step oko nominalnog
% 
x10 = Xe;
x20 = Se;

delta = S_upper_limit-Se;
S_step1 = delta;
S_step2 = delta;
S_time_step1 = 20;
S_time_step2 = 50;

disp(['Ravnotežna koncentracija supstrata: ',num2str(Se)])
disp(['Inicijalna koncentracija supstrata: ', num2str(x20)])
disp(['Step koncentracije susptrata na ulazu: ', num2str(S_step1)])

%% Projektovanje linearnog kontrolera 
w0 = 0.5637; % rad/h, očitano sa grafika

w1_design = w0;
s = tf('s');
%% PI regulator
Fpf_design = 90;
G_design = G;
Ti1 = 1/w1_design*tan(Fpf_design-pi/2-unwrap(angle(freqresp(G_design,w1_design))));
Kc1 = 1/abs(freqresp((1+1/Ti1/s)*G_design,w1_design));

KPI = Kc1*(1 + 1/Ti1/s);

KPI_aw = 1/(Ti1*s+1);
Tf = 1/(20*w1_design);

Td1 = Ti1/2; % eksperimentalno procenjeno
Kfb =0;%Kc1*Td1*s/(Tf*s + 1);

%% Uvođenje šuma merenja
std_noise = 0;%1e-5*Se;
disp(['šum pri merenju: ',num2str(std_noise)]);
%% Uvođenje poremećaja
disturb_X = 0;%0.2*Xe;
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);
%% Pokretanje simulacije
sim(['robusni_lin_regulacija_nelin_sim.slx']);
%% Rezultati
%PROMENLJIVE STANJA U VREMENU
figure;
subplot(2,1,1)
hold all;
plot(t_out,x1_out)
plot(t_out,x1_out1)
plot(t_out,x1_out2)
%plot(t_out,x1_out3)
%plot(t_out,x1_out4)
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
%legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
grid

subplot(2,1,2)
hold all;
plot(t_out,x2_out)
plot(t_out,x2_out1)
plot(t_out,x2_out2)
%plot(t_out,x2_out3)
%plot(t_out,x2_out4)
xlabel('vreme [h]')
ylabel('x2 = S [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
%legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
sgtitle('Promenljive stanja')
grid

%% 

%FAZNI DIJAGRAM

figure;
hold all;
plot(x1_out,x2_out)
plot(x1_out1,x2_out1)
plot(x1_out2,x2_out2)
%plot(x1_out3,x2_out3)
%plot(x1_out4,x2_out4)
xlabel('x1 = X [g/l]')
ylabel('x2 = S [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
%legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
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
%plot(t_out,F_out3)
%plot(t_out,F_out4)
grid
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
%legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
title('Zapremisnki protok kroz reaktor')
ylabel('u = F [l/h]')
xlabel('vreme [h]')

subplot(2,1,2)
hold all;
plot(t_out,S_out)
plot(t_out,S_out1)
plot(t_out,S_out2)
%plot(t_out,S_out3)
%plot(t_out,S_out4)
grid
ylabel('y = S [g/l]')
xlabel('vreme [h]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
%legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}','1.1Y_{nom}','1.2Y_{nom}')
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
xlabel('vreme [h]')
ylabel('x1 = X [g/l]')
legend('Y_{nom}','0.8Y_{nom}','0.9Y_{nom}')
title('Koncentracija biomase na izlazu')
grid

