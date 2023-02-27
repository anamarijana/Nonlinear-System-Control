clc
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


%% Limit upravljanja 
F_upper_limit = Fe+1e-3*Fe;
F_lower_limit = 0;

% procenjeno sa grafika
S_upper_limit = 0.2325;
S_lowe_limit = 0;
%% Ograničenje propusnog opsega
w0 = 0.5637; % rad/h, očitano sa grafika

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
S_time_step2 = 40;
%% Uvođenje šuma merenja
std_noise =0;%1e-5*Se; %chattering
%% Uvođenje poremećaja
disturb_X = 0.2*Xe;
disturb_X_time = 5;
%% Klizno upravljanje proporcionalno
beta = S_upper_limit*w0/2;
Ki = 0;
%% uvođenje graničnog sloja

% beta/Fi ne treba da bude dominantno beta/Fi = 10*w0
Fi = beta/5/w0; 
%% Pokretanje simulacije
sim(['SMC_regulacija_nelin_sim.slx']);
%% Procena propusnog opsega tranzijent
initial = x20 - Se;
if(initial)
    response = -S_out+S_out(1);
    yend = response(end);
    idx = find(response<=0.63*abs(yend));
    max_idx = idx(end);
    time_constant = t_out(max_idx);
    w0_estimated = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated)
    figure;
    grid
    movegui('center')
    hold all;
    plot(t_out,response)
    plot(t_out(max_idx),response(max_idx),'r.')
    title(['Tranzijent, procenjen propusni opseg: ' num2str(w0_estimated)])
end
    
%% Procena propusnog opsega step
if(S_step)
    steady_state_time = S_time_step2-5;
    steady_state_index = find(t_out==steady_state_time);
    response1 = S_out(1:steady_state_index);
    t_response = t_out(1:steady_state_index);
    response = response1-response1(1);
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_out(max_idx)-S_time_step1;
    w0_estimated = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated)
    figure;
    grid
    movegui('center')
    hold all;
    plot(t_response,response)
    plot(t_response(max_idx),response(max_idx),'r.')
    title(['Step: procenjen propusni opseg; ' num2str(w0_estimated)])
end

%% Procena propusnog opsega poremecaj
if(disturb_X)
    [value_start,index_start] = min(S_out);
    response1 = S_out(index_start:end);
    t_response = t_out(index_start:end);
    response  = response1 - value_start;
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_response(max_idx)-disturb_X_time;
    w0_estimated = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated)
    figure;
    grid
    movegui('center')
    hold all;
    plot(t_response,response)
    plot(t_response(max_idx),response(max_idx),'r.')
    title(['Poremecaj: procenjen propusni opseg; ' num2str(w0_estimated)])
end
%% Preskok izlaza
preskok = (max(S_out)-S_out(end))/S_out(end);
disp(['Preskok koncentracije supstrata: ',num2str(preskok)]);

%% Greška stacionarnog stanja izlaza 
e_stac_stanje = S_out(end)-Se;
disp(['Greška stacionarnog stanja izlaza: ',num2str(e_stac_stanje)]);
%% Vremenske zavisnosti
figure;
movegui('southwest')
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
movegui('south')
plot(x1_out,x2_out)
xlabel('x1 [g/l]')
ylabel('x2 [g/l]')
grid
title('regulacija: fazni dijagram')

%%
figure;
movegui('southeast')
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

%% 
disp(['Ravnotežna koncentracija supstrata: ',num2str(Se)])
disp(['Inicijalna koncentracija supstrata: ', num2str(x20)])
disp(['Step koncentracije susptrata na ulazu: ', num2str(S_step)])
disp(['šum pri merenju: ',num2str(std_noise)]);
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);