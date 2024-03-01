clc;
clear;
close all;


avg_N = 1000;                           % Average times
variance_arr = [0.08 0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = [round(power(10,0),5) round(power(10,0.5),5) round(power(10,1),5) round(power(10,1.5),5) round(power(10,2),5)];
% Figure 2 with different transmit power

avg_mse_nonrobust = zeros(4,5);
avg_ser_nonrobust = zeros(4,5);
iteration = 1000;                       %parfor iteration

%% Simulation start

%Progressbar
ppm = ParforProgressbar(iteration);
ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


parfor i=1:1000

    % Data for figure2
    temp_nonrobust_ser = zeros(4,5);
    temp_nonrobust = zeros(4,5);
    for sig = 1:length(variance_arr)
        for power_index = 1:length(tx_power_arr)
            [mse,ser,power_array]=mmse(40,variance_arr(sig),tx_power_arr(power_index),2,0);
            temp_nonrobust(sig,power_index) = mse;
            temp_nonrobust_ser(sig,power_index) = ser;
        end
    end
    avg_mse_nonrobust = avg_mse_nonrobust+temp_nonrobust;
    avg_ser_nonrobust = avg_ser_nonrobust+temp_nonrobust_ser;



    %Progressbar
    pause(100/iteration);
    ppm.increment();
end

%Delete Progressbar
delete(ppm);

%% Plot for figure2

%Average result
avg_mse_nonrobust = avg_mse_nonrobust./avg_N;
avg_ser_nonrobust = avg_ser_nonrobust./avg_N;

%Change into dBm
tx_power_arr = 10*log10(tx_power_arr);

array1 = [0.0107    0.0102    0.0079    0.0074    0.0068];  % Robust(paper), variance=0.05
array2 = [0.0059    0.0045    0.0032    0.0025    0.0022];  % Robust(paper), variance=0.01
array3 = [0.0050    0.0035    0.0022    0.0020    0.0017];  % Robust(paper), variance=0.002

figure(2)
semilogy(tx_power_arr,avg_mse_nonrobust(1,:),'-o',tx_power_arr,avg_mse_nonrobust(2,:),'-square',tx_power_arr,avg_mse_nonrobust(3,:),'-*',tx_power_arr,avg_mse_nonrobust(4,:),'-^',tx_power_arr,array1,'-',tx_power_arr,array2,'-',tx_power_arr,array3,'-');
colororder([1 0 0;0 0 1;0 0 0;1 0 1; 0 0 1;0 0 0;1 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.08','Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002','Robust(paper), {\sigma^2}=0.05','Robust(paper), {\sigma^2}=0.01','Robust(paper), {\sigma^2}=0.002')



noise_var = power(10,-11);
N0 = noise_var;
M = 16;
k=log2(M);
for i=1:length(tx_power_arr)
    SNR=exp(tx_power_arr(i)*log(10)/10);    	% signal-to-noise ratio
    % theoretical symbol error rate
    theo_err_prb(i)=4*qfunc(sqrt(3*k*SNR/(M-1)));
end

figure(3)
semilogy(tx_power_arr,avg_ser_nonrobust(1,:),'-o',tx_power_arr,avg_ser_nonrobust(2,:),'-square',tx_power_arr,avg_ser_nonrobust(3,:),'-*',tx_power_arr,theo_err_prb,'--');
colororder([1 0 0;0 0 1;0 0 0;1 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("SER")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002','Theoretical')
