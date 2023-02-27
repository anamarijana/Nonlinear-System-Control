clc;
clear;
close all;
%% Parametri modela

tfin=50;     % sec

V = 4; %l
SF = 10; %g/l
Y = 0.5;
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
G_lin_poles = pole(G);
G_lin_zeros = zero(G);

%% Ponašanje oko nominalnog režima

S_initial = Se+0.1*Se; %proizvoljna tačka dovoljno blizu
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);

% x10 = X_initial;
% x20 = S_initial;
% 
x10 = Xe;
x20 = Se;

S_step = 0;%0.01*Se;
S_time_step1 = 2;
S_time_step2 = 35;


%% Limit upravljanja 

F_upper_limit = 1e-3*Fe;
F_lower_limit = 0;
%% Projektovanje linearnog kontrolera 
w0 = 0.5637; % rad/h
w1_design = 0.5637;
s = tf('s');

Fpf_design = 90;
G_design = G;
Ti1 = 1/w1_design*tan(Fpf_design-pi/2-unwrap(angle(freqresp(G_design,w1_design))));
Kc1 = 1/abs(freqresp((1+1/Ti1/s)*G_design,w1_design));
KPI = Kc1*(1 + 1/Ti1/s);

KPI_aw = 1/(Ti1*s+1);
Tf = 1/(20*w1_design);

Td1 = Ti1/2;
Kfb1 = Kc1*Td1*s/(Tf*s + 1);

Td2 = Ti1/3;
Kfb2 = Kc1*Td2*s/(Tf*s + 1);

Td3 = Ti1/4;
Kfb3 = Kc1*Td3*s/(Tf*s + 1);

Td4 = Ti1/5;
Kfb4 = Kc1*Td4*s/(Tf*s + 1);

Td5 = Ti1/6;
Kfb5 = Kc1*Td5*s/(Tf*s + 1);

%% Uvođenje šuma merenja
std_noise = 1e-5*Se;

%% Uvođenje poremećaja

disturb_X = 0.2*Xe;
%% Pokretanje simulacije
sim('diferencijalni_lin_regulacija_nelin_sim.slx');

%% Vremenske zavisnosti

figure;
hold all;
plot(t_out,x1_out)
plot(t_out,x1_out1)
plot(t_out,x1_out2)
plot(t_out,x1_out3)
plot(t_out,x1_out4)
plot(t_out,x1_out5)
xlabel('vreme [h]')
ylabel('X [g/l]')
grid
legend('PI regulator','Td = Ti/2','Td = Ti/3','Td = Ti/4','Td = Ti/5','Td = Ti/6')

figure;
hold all;
plot(t_out,x2_out)
plot(t_out,x2_out1)
plot(t_out,x2_out2)
plot(t_out,x2_out3)
plot(t_out,x2_out4)
plot(t_out,x2_out5)
xlabel('vreme [h]')
ylabel('S [g/l]')
grid
legend('PI regulator','Td = Ti/2','Td = Ti/3','Td = Ti/4','Td = Ti/5','Td = Ti/6')
sgtitle('Promenljive stanja u vremenu')

%%
figure;
hold all;
plot(t_out,F_out)
plot(t_out,F_out1)
% plot(t_out,F_out2)
% plot(t_out,F_out3)
% plot(t_out,F_out4)
plot(t_out,F_out5)
grid
title('Protok supstrata u reaktor')
ylabel('F [l/h]')
xlabel('vreme [h]')
legend('PI regulator','Td = Ti/2','Td = Ti/6')

%legend('PI regulator','Td = Ti/2','Td = Ti/3','Td = Ti/4','Td = Ti/5','Td = Ti/6')

%%
figure;
hold all;
plot(t_out,S_out)
plot(t_out,S_out1)
plot(t_out,S_out5)
plot(t_out,S_ref_out,'k--')
grid
ylabel('S [g/l]')
xlabel('vreme [h]')
legend('PI regulator','Td = Ti/2','Td = Ti/6')
%legend('PI regulator','Td = Ti/2','Td = Ti/3','Td = Ti/4','Td = Ti/5','Td = Ti/6')
title('Koncentracija supstrata na izlazu')


