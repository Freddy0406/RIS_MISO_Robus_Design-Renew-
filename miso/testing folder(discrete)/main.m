clc;
clear;
close all;


avg_N = 1000;                           % Average times
variance_arr = [0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 2 with different transmit power

avg_mse_robust = zeros(3,5);
avg_ser_robust = zeros(3,5);
power_list = zeros(5,1000,10000);
iteration = 1000;                       %parfor iteration

%% Simulation start

%Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


parfor i=1:1000

    % Data for figure2
    temp_robust = zeros(3,5);
    temp_robust_ser = zeros(3,5);
    temp_nonrobust = zeros(3,5);
    temp_power = zeros(5,1000,10000);
    for sig = 1:length(variance_arr)
        for power_index = 1:length(tx_power_arr)
            [mse,ser,power_array]=mmse(40,variance_arr(sig),tx_power_arr(power_index),3,2);
            temp_robust(sig,power_index) = mse;
            temp_robust_ser(sig,power_index) = ser;
            temp_power(power_index,i,:) = power_array;
        end
    end
    avg_mse_robust = avg_mse_robust+temp_robust;
    avg_ser_robust = avg_ser_robust+temp_robust_ser;
    power_list = power_list+temp_power;



    %Progressbar
    pause(100/iteration);
    ppm.increment();
end

%Delete Progressbar
delete(ppm);

%% Plot for figure2

%Average result
avg_mse_robust = avg_mse_robust./avg_N;
avg_ser_robust = avg_ser_robust./avg_N;

%Change into dBm
tx_power_arr = 10*log10(tx_power_arr);

figure(2)
semilogy(tx_power_arr,avg_mse_robust(1,:),'-o',tx_power_arr,avg_mse_robust(2,:),'-square',tx_power_arr,avg_mse_robust(3,:),'-*');
colororder([1 0 0;0 0 1;0 0 0]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002')


figure(3)
semilogy(tx_power_arr,avg_ser_robust(1,:),'-o',tx_power_arr,avg_ser_robust(2,:),'-square',tx_power_arr,avg_ser_robust(3,:),'-*');
colororder([1 0 0;0 0 1;0 0 0]);
xlabel("Transmit Power (dBm)")
ylabel("SER")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002')
