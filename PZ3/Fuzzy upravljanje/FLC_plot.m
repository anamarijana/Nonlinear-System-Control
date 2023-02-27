clc;
clear;
close all;
%%
FLC= readfis('last');

%% Projektovani FLC: Pravila
showrule(FLC) 
%% MFS input 1
input1_mfs = FLC.Inputs(1).MembershipFunctions;
num_of_mfs = length(input1_mfs);
pulse = [0:1, 0:-1:0];
time1 = [-3 -2 0];
time2 = [-0.5 0 0.5];
time3 = [0 2 3];
figure; 
hold all;
plot(time1, pulse)
plot(time2, pulse)
plot(time3, pulse)
pera = [-1.9 0 1.9];
for i=1:num_of_mfs
    text(pera(i),1.1,input1_mfs(i).Name);
end
xticks([-3 -2 -0.5 0 0.5 2 3])
xticklabels({'-0.1299','-0.0137','-0.00025','0','0.00025','0.2187','0.3349'})
xlim([-2 2])
title('Fuzzy skupovi za ulaz e');
xlabel('e'); ylabel('Stepen pripadnosti'); a=axis; a(4)=1.2; axis(a);
% plotmf(FLC,'input',1); title('Input 1 = Greska e')
%% MFS input 2

input2_mfs = FLC.Inputs(2).MembershipFunctions;
num_of_mfs = length(input2_mfs);
pulse = [0:1, 0:-1:0];
time1 = [-14.4203 -7.8663 -1.3123];
time2 = [-6.554 0 6.554];
time3 = [1.3123 7.8663 14.4203];
figure; 
hold all;
plot(time1, pulse)
plot(time2, pulse)
plot(time3, pulse)
pera = [-7.5 0 7.5];
for i=1:num_of_mfs
    text(pera(i),1.1,input2_mfs(i).Name);
end
xticks([-14.42 -7.87 -6.55 -1.31 0 1.31 6.55 7.87 14.42])
%xticklabels({'-14.4203','-0.0137','-0.00025','0','0.00025','0.2187','0.3349'})
xlim([-7.8663 7.8663])
title('Fuzzy skupovi za ulaz ed');
xlabel('ed'); ylabel('Stepen pripadnosti'); a=axis; a(4)=1.2; axis(a);
% plotmf(FLC,'input',1); title('Input 1 = Greska e')
%% MFS output
output_singleton_mfs = FLC.Output.MembershipFunctions;
num_of_singletons = length(output_singleton_mfs);
figure;
hold all;
for i=1:num_of_singletons
    stem(i-1, 1);
    text(i-1.2,1.1,output_singleton_mfs(i).Name);
end
xlim([-1 5])
xticks([0 1 2 3 4])
xticklabels({'0','1.6044','3.2089','3.2105','3.2121'})
title('Fuzzy skupovi za izlaz u');
xlabel('u'); ylabel('Stepen pripadnosti'); a=axis; a(4)=1.2; axis(a);
%% Projektovani FLC: Ulazno-izlazna povrsina
figure; gensurf(FLC);
