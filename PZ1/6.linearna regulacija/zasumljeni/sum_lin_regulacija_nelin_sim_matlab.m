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
x10 = X_initial;
x20 = S_initial;


%% step oko nominalnog
% 
% x10 = Xe;
% x20 = Se;

delta = 0;%0.2325-Se;
S_step1 = delta;
S_step2 = delta;
S_time_step1 = 5;
S_time_step2 = 40;

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
std_noise = 1e-5*Se;
disp(['šum pri merenju: ',num2str(std_noise)]);
%% Uvođenje poremećaja
disturb_X = 0.2*Xe;
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);
%% Pokretanje simulacije
sim(['sum_lin_regulacija_nelin_sim.slx']);
%% Procena propusnog opsega
response = -S_out+S_out(1);
yend = response(end);
idx = find(response<=0.63*abs(yend));
max_idx = idx(end);
time_constant = t_out(max_idx);
w0_estimated = 1./time_constant;
disp('Propusni opseg regulacije:')
disp(w0_estimated)
%% Vremenske zavisnosti
figure;
subplot(2,1,1)
plot(t_out,x1_out)
xlabel('vreme [h]')
ylabel('X [g/l]')
grid

subplot(2,1,2)
plot(t_out,x2_out)
xlabel('vreme [h]')
ylabel('S [g/l]')
grid
sgtitle('PI regulacija: promenljive stanja u vremenu')

%% 

figure;
plot(x1_out,x2_out)
xlabel('x1 [g/l]')
ylabel('x2 [g/l]')
grid
title('PI regulacija: fazni dijagram')

%%
figure;
subplot(2,1,1)
hold all;
plot(t_out,F_out)
grid
title('Protok supstrata u reaktor')
ylabel('F [l/h]')
xlabel('vreme [h]')

subplot(2,1,2)
hold all;
plot(t_out,S_out)
plot(t_out,S_ref_out,'k--')
grid
ylabel('S [g/l]')
xlabel('vreme [h]')
legend('PI regulator','referentna vrednost')
title('Koncentracija supstrata na izlazu')


