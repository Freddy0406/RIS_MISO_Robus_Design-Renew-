clc;
clear;
close all;


avg_N = 1000;                           % Average times
N_arr = [10 20 30 40 50 60];
% Figure 2 with different transmit power

avg_mse_robust = zeros(3,6);
avg_ser_robust = zeros(3,6);
power_list = zeros(6,1000,10000);
iteration = 1000;                       %parfor iteration

%% Simulation start

%Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


parfor i=1:1000

    % Data for figure2
    temp_robust = zeros(3,6);
    temp_robust_ser = zeros(3,6);
    temp_power = zeros(6,1000,10000);
    for N_index = 1:length(N_arr)
        [mse,ser,power_array]=mmse(N_arr(N_index),0,round(power(10,1),5),1,0);
        temp_robust(1,N_index) = mse;
        temp_robust_ser(1,N_index) = ser;
        temp_power(N_index,i,:) = power_array;
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
figure(2)
semilogy(N_arr,avg_mse_robust(1,:),'-o');
colororder([1 0 0]);
xlabel("Reflecting elements")
ylabel("MSE")
legend('Perfect')


% figure(3)
% semilogy(tx_power_arr,avg_ser_robust(1,:),'-o',tx_power_arr,avg_ser_robust(2,:),'-square',tx_power_arr,avg_ser_robust(3,:),'-*');
% colororder([1 0 0;0 0 1;0 0 0]);
% xlabel("Transmit Power (dBm)")
% ylabel("SER")
% legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002')
