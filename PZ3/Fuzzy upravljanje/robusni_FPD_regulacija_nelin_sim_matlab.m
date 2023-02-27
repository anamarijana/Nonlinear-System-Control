
clear;
close all;
%%
sim_file_name = 'robusni_FPD_regulacija_nelin_sim';
FLC= readfis('last'); 
%% Parametri modela

tfin=60;     % sec
V = 4.0000; %l
SF = 10.0000; %g/l
Y = 0.5000;
mi_max = 1;
K1 = 0.0300; %g/l
K2 = 0.5000; %l/g

w0 = 0.5637;

%% Ispitivanje robusnosti
Y1 = Y - 0.2*Y;
Y2 = Y - 0.1*Y;
Y3 = Y + 0.1*Y;
Y4 = Y + 0.2*Y;

%% Nominalni režim, ranije izračunato

Se = 0.218662692463450; 
Xe = 4.89066865376828;
Fe = 3.20891060125344;
mie = mi_max*Se/(K2*Se^2 + Se + K1);
%% Limiti
F_upper_limit = Fe + 1e-3*Fe;
F_lower_limit = 0;

S_upper_limit = 0.2324;
S_lower_limit = 0;

derivative_S_upper_limit = 7.866;
derivative_S_lower_limit = -7.866;

e_upper_limit = 0.2187;
e_lower_limit = -0.0137;

derivative_e_upper_limit = 7.8663;
derivative_e_lower_limit = -7.8663;

derivative_F_upper_limit = 3.2121e+04;
derivative_F_lower_limit = -3.2121e+04;

%% Uklanjanje početnog uslova

S_initial = S_upper_limit; %proizvoljna tačka dovoljno blizu
mi_initial = mi_max*S_initial/(K2*S_initial^2 + S_initial + K1); 
F_initial = mi_initial*V;
X_initial = Y*(SF-S_initial);
% 
% x10 = X_initial;
% x20 = S_initial;
% mi0 = mi_initial;

%% Stepovi oko nominalnog režima
x10 = Xe;
x20 = Se;
mi0 = mie;

delta = 0;% S_upper_limit-Se;
S_step1 = delta;
S_step2 = delta;
S_time_step1 = 20;
S_time_step2 = 50;

%% Uvođenje šuma merenja
std_noise = 0;%1e-5*Se;
disp(['šum pri merenju: ',num2str(std_noise)]);
%% Uvođenje poremećaja
disturb_X =  0.2*Xe;
disturb_X_time = 5;
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);

%% Pokretanje simulacije
open(sim_file_name)
sim(sim_file_name);
%% Rezultati
%PROMENLJIVE STANJA U VREMENU
figure;
movegui('southwest')
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
movegui('south')
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
movegui('southeast')
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

%% 
disp(['Ravnotežna koncentracija supstrata: ',num2str(Se)])
disp(['Inicijalna koncentracija supstrata: ', num2str(x20)])
disp(['Step koncentracije susptrata na ulazu: ', num2str(S_step1)])
disp(['Poremećaj koncentracije biomase: ',num2str(disturb_X)]);
disp(['šum pri merenju: ',num2str(std_noise)]);


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
    
    response = S_out1;
    yend = response(end);
    idx = find(response<=0.63*abs(yend));
    max_idx = idx(end);
    time_constant = t_out(max_idx);
    w0_estimated1 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated1)
    figure;
    grid
    movegui('center')
    hold all;
    plot(t_out,response)
    plot(t_out(max_idx),response(max_idx),'r.')
    title(['Tranzijent, procenjen propusni opseg: ' num2str(w0_estimated1)])
    
    response = S_out2;
    yend = response(end);
    idx = find(response<=0.63*abs(yend));
    max_idx = idx(end);
    time_constant = t_out(max_idx);
    w0_estimated2 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated2)
    figure;
    grid
    movegui('center')
    hold all;
    plot(t_out,response)
    plot(t_out(max_idx),response(max_idx),'r.')
    title(['Tranzijent, procenjen propusni opseg: ' num2str(w0_estimated2)])
    
    response = -S_out3+S_out3(1);
    yend = response(end);
    idx = find(response<=0.63*abs(yend));
    max_idx = idx(end);
    time_constant = t_out(max_idx);
    w0_estimated3 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated3)
    
    response = -S_out4+S_out4(1);
    yend = response(end);
    idx = find(response<=0.63*abs(yend));
    max_idx = idx(end);
    time_constant = t_out(max_idx);
    w0_estimated4 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated4)
    
end


%% Procena propusnog opsega step
if(S_step1)
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
    
    steady_state_time = S_time_step2-5;
    steady_state_index = find(t_out==steady_state_time);
    response1 = S_out1(1:steady_state_index);
    t_response = t_out(1:steady_state_index);
    response = response1-response1(1);
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_out(max_idx)-S_time_step1;
    w0_estimated1 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated1)
    
    steady_state_time = S_time_step2-5;
    steady_state_index = find(t_out==steady_state_time);
    response1 = S_out2(1:steady_state_index);
    t_response = t_out(1:steady_state_index);
    response = response1-response1(1);
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_out(max_idx)-S_time_step1;
    w0_estimated2 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated2)
    
    steady_state_time = S_time_step2-5;
    steady_state_index = find(t_out==steady_state_time);
    response1 = S_out3(1:steady_state_index);
    t_response = t_out(1:steady_state_index);
    response = response1-response1(1);
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_out(max_idx)-S_time_step1;
    w0_estimated3 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated3)
    
    steady_state_time = S_time_step2-5;
    steady_state_index = find(t_out==steady_state_time);
    response1 = S_out4(1:steady_state_index);
    t_response = t_out(1:steady_state_index);
    response = response1-response1(1);
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_out(max_idx)-S_time_step1;
    w0_estimated4 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated4)

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
    
    [value_start,index_start] = min(S_out);
    response1 = S_out1(index_start:end);
    t_response = t_out(index_start:end);
    response  = response1 - value_start;
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_response(max_idx)-disturb_X_time;
    w0_estimated1 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated1)
    
    [value_start,index_start] = min(S_out);
    response1 = S_out2(index_start:end);
    t_response = t_out(index_start:end);
    response  = response1 - value_start;
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_response(max_idx)-disturb_X_time;
    w0_estimated2 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated2)
    
    [value_start,index_start] = min(S_out);
    response1 = S_out3(index_start:end);
    t_response = t_out(index_start:end);
    response  = response1 - value_start;
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_response(max_idx)-disturb_X_time;
    w0_estimated3 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated3)
    
    [value_start,index_start] = min(S_out);
    response1 = S_out4(index_start:end);
    t_response = t_out(index_start:end);
    response  = response1 - value_start;
    yend = response(end);
    idx = find(response<=0.63*yend);
    max_idx = idx(end);
    time_constant = t_response(max_idx)-disturb_X_time;
    w0_estimated4 = 1./time_constant;
    disp('Propusni opseg regulacije:')
    disp(w0_estimated4)
     
end