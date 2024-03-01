clc;
clear;
close all;

%test
avg_N = 1000;                           % Average times
variance_arr = [0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = zeros(1,20);
%Figure 2 with different transmit power
%Let tx_power_arr = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20] in dBm
%Change to linear
dBm_index = 0;

% for z = 1:11
%     tx_power_arr(z) = round((power(10,dBm_index/10)),5);
%     dBm_index = dBm_index+5;
% end
% dBm_index = 12;
% 
% for z = 12:length(tx_power_arr)
%     tx_power_arr(z) = round((power(10,dBm_index/10)),5);
%     dBm_index = dBm_index+10;
% end

for z = 1:length(tx_power_arr)
    tx_power_arr(z) = round((power(10,dBm_index/10)),5);
    dBm_index = dBm_index+2;
end

avg_mse_robust = zeros(1,length(tx_power_arr));
avg_mse_robust_wn = zeros(1,length(tx_power_arr));
avg_ser_robust = zeros(1,length(tx_power_arr));
power_list = zeros(length(tx_power_arr),1000,10000);
iteration = 1000;                       %parfor iteration

%% Simulation start

%Progressbar
% ppm = ParforProgressbar(iteration);
% ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


for i=1:1000

    % Data for figure2
    temp_robust = zeros(1,length(tx_power_arr));
    temp_robust_wn = zeros(1,length(tx_power_arr));
    temp_robust_ser = zeros(1,length(tx_power_arr));
    temp_power = zeros(length(tx_power_arr),1000,10000);
    for power_index = 1:length(tx_power_arr)
        [mse,ser,power_array]=mmse(40,0,tx_power_arr(power_index),1,0);
        temp_robust(1,power_index) = mse;
        temp_robust_ser(1,power_index) = ser;
        [mse,~,~]=mmse(40,0,tx_power_arr(power_index),2,0);
        temp_robust_wn(1,power_index) = mse;
        temp_power(power_index,i,:) = power_array;
    end

    avg_mse_robust = avg_mse_robust+temp_robust;
    avg_mse_robust_wn = avg_mse_robust_wn+temp_robust_wn;
    avg_ser_robust = avg_ser_robust+temp_robust_ser;
    power_list = power_list+temp_power;


    fprintf('\n\t%d time completed...\n\n',i);
    %Progressbar
    % pause(100/iteration);
    % ppm.increment();
end

%Delete Progressbar
% delete(ppm);

%% Plot for figure2

%Average result
avg_mse_robust = avg_mse_robust./avg_N;
avg_mse_robust_wn = avg_mse_robust_wn./avg_N;
avg_ser_robust = avg_ser_robust./avg_N;

%Change into dBm
tx_power_arr = 10*log10(tx_power_arr);

figure(2)
semilogy(tx_power_arr,avg_mse_robust(1,:),'-o',tx_power_arr,avg_mse_robust_wn(1,:),'-o');
colororder([1 0 0;0 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Perfect channel','Perfect channel(without norm)')

figure(2)
semilogy(tx_power_arr,avg_mse_robust_wn(1,:),'-o');
colororder([0 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Perfect channel')


noise_var = power(10,-11);
N0 = noise_var;
M = 16;
k=log2(M);
theo_err_prb = zeros(1,length(tx_power_arr));
for i=1:length(tx_power_arr)
    SNR=exp(tx_power_arr(i)*log(10)/10);    	% signal-to-noise ratio
    % theoretical symbol error rate
    theo_err_prb(i)=4*qfunc(sqrt(3*k*SNR/(M-1)));
end

figure(3)
semilogy(tx_power_arr,avg_ser_robust(1,:),'-o',tx_power_arr,theo_err_prb,'--');
colororder([1 0 0;1 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("SER")
legend('Perfect channel','Theoretical')
