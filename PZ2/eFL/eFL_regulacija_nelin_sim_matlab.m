clc;
clear;
close all;
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


% procenjeno sa grafika
S_upper_limit = 0.2325;
S_lowe_limit = 0;

%% Otklanjanje početnog uslova

S_initial = S_upper_limit;
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);
% 
% x10 = X_initial;
% x20 = S_initial;
%% Step odziv reference
% 
x10 = Xe;
x20 = Se;

S_step = 0;%S_upper_limit-Se;
S_time_step1 = 2;
S_time_step2 = 40; %80 racunanje w0 da se ne menja formula

%% 
disp(['Ravnotežna koncentracija supstrata: ',num2str(Se)])
disp(['Inicijalna koncentracija supstrata: ', num2str(x20)])
disp(['Step koncentracije susptrata na ulazu: ', num2str(S_step)])

%% Limit upravljanja 
F_upper_limit = Fe+1e-3*Fe;
F_lower_limit = 0;
%% Ograničenje propusnog opsega
w0 = 0.5637; % rad/h, očitano sa grafika
%% Uvođenje šuma merenja
std_noise = 1e-5*Se;
disp(['šum pri merenju: ',num2str(std_noise)]);
%% Uvođenje poremećaja
disturb_X = 0.2*Xe;
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);

disturb_S = 0;%Se;
disp(['Poremećaj koncentracije supstrata: ',num2str(disturb_S)]);

%% Linearizacija
K0_FLRT = w0; %FLRT
syms x
f_fli_rt_zelj = collect((x+2*w0)^2, x);

K0_FLIRT = 2.2548;
Ki = 1.27;
%% Pokretanje simulacije
open('eFL_regulacija_nelin_sim.slx')
sim(['eFL_regulacija_nelin_sim.slx']);
%% Procena propusnog opsega 
response = -S_out+S_out(1);
yend = response(end); %tranzijent i stacionarno
idx = find(response<=0.63*abs(yend));
max_idx = idx(end);
time_constant = t_out(max_idx);
w0_estimated = 1./time_constant;
disp('Propusni opseg regulacije:')
disp(w0_estimated)
%% Preskok izlaza

preskok = (max(S_out)-S_out(end))/S_out(end);
disp(['Preskok koncentracije supstrata: ',num2str(preskok)]);

%% Greška stacionarnog stanja izlaza 
e_stac_stanje = S_out(end)-Se;
disp(['Greška stacionarnog stanja izlaza: ',num2str(e_stac_stanje)]);

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
sgtitle('regulacija: promenljive stanja u vremenu')

%% 

figure;
plot(x1_out,x2_out)
xlabel('x1 [g/l]')
ylabel('x2 [g/l]')
grid
title('regulacija: fazni dijagram')

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
legend('regulator','referentna vrednost')
title('Koncentracija supstrata na izlazu')


