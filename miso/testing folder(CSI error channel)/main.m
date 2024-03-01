clc;
clear;
close all;



variance_arr = [0.05 0.01 0.002];       % Figure 1 with different variance
tx_power_arr = zeros(1,50);
%Figure 2 with different transmit power
%Let tx_power_arr = [0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20] in dBm
%Change to linear
dBm_index = 0;

% for z = 1:11
%     tx_power_arr(z) = round((power(10,dBm_index/10)),5);
%     dBm_index = dBm_index+1;
% end
% dBm_index = 12;
% 
% for z = 12:length(tx_power_arr)
%     tx_power_arr(z) = round((power(10,dBm_index/10)),5);
%     dBm_index = dBm_index+2;
% end

for z = 1:length(tx_power_arr)
    tx_power_arr(z) = round((power(10,dBm_index/10)),5);
    dBm_index = dBm_index+2;
end

avg_mse_robust = zeros(3,length(tx_power_arr));
avg_mse_robust_wn = zeros(3,length(tx_power_arr));
avg_ser_robust = zeros(3,length(tx_power_arr));
avg_N = 1000;                           % Average times
iteration = 1000;                       %parfor iteration

%% Simulation start

%Progressbar
% ppm = ParforProgressbar(iteration);
% ppm = ParforProgressbar(iteration, 'showWorkerProgress', true);


for i=1:1000
    % Data for figure2
    temp_robust = zeros(3,length(tx_power_arr));
    temp_robust_wn = zeros(3,length(tx_power_arr));
    temp_robust_ser = zeros(3,length(tx_power_arr));
    temp_nonrobust = zeros(3,length(tx_power_arr));
    for sig = 1:length(variance_arr)
        for power_index = 1:length(tx_power_arr)
%             [mse,ser,power_array]=mmse(40,variance_arr(sig),tx_power_arr(power_index),1,0);
%             temp_robust(sig,power_index) = mse;
%             temp_robust_ser(sig,power_index) = ser;
            [mse,~,~]=mmse(100,variance_arr(sig),tx_power_arr(power_index),2,0);
            temp_robust_wn(sig,power_index) = mse;
        end
    end
    avg_mse_robust = avg_mse_robust+temp_robust;
    avg_mse_robust_wn = avg_mse_robust_wn+temp_robust_wn;
    avg_ser_robust = avg_ser_robust+temp_robust_ser;

    fprintf('\t%d time completed\n',i);

    %Progressbar
%     pause(100/iteration);
%     ppm.increment();
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

array1 = [10^-1.91585    10^-1.98536    10^-2.05    10^-2.06429    10^-2.085714];  % Robust(paper), variance=0.05
array2 = [10^-2.1125    10^-2.18375    10^-2.3    10^-2.405    10^-2.465];  % Robust(paper), variance=0.01
array3 = [10^-2.15    10^-2.26357     10^-2.47    10^-2.72812    10^-2.92875];  % Robust(paper), variance=0.002

tx_power_arr_paper = [0,5,10,15,20];

figure(2)
semilogy(tx_power_arr,avg_mse_robust(1,:),'-o',tx_power_arr,avg_mse_robust(2,:),'-square',tx_power_arr,avg_mse_robust(3,:),'-d',tx_power_arr_paper,array1,'-*',tx_power_arr_paper,array2,'-',tx_power_arr_paper,array3,'-',tx_power_arr,avg_mse_robust_wn(1,:),'--o',tx_power_arr,avg_mse_robust_wn(2,:),'--square',tx_power_arr,avg_mse_robust_wn(3,:),'--d');
grid on;
% xlim([0 20])
colororder([1 0 0;0 0 1;0 0 0; 1 0 0;0 0 1;0 0 0; 1 0 0;0 0 1;0 0 0]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002','Robust(paper), {\sigma^2}=0.05','Robust(paper), {\sigma^2}=0.01','Robust(paper), {\sigma^2}=0.002','Robust(without norm), {\sigma^2}=0.05','Robust(without norm), {\sigma^2}=0.01','Robust(without norm), {\sigma^2}=0.002')


figure(2)
semilogy(tx_power_arr_paper,array1,'-',tx_power_arr_paper,array2,'-',tx_power_arr_paper,array3,'-',tx_power_arr,avg_mse_robust_wn(1,:),'--o',tx_power_arr,avg_mse_robust_wn(2,:),'--square',tx_power_arr,avg_mse_robust_wn(3,:),'--d');
grid on;
xlim([0 20])
colororder([1 0 0;0 0 1;0 0 0; 1 0 0;0 0 1;0 0 0]);
xlabel("Transmit Power (dBm)")
ylabel("MSE")
legend('Robust(paper), {\sigma^2}=0.05','Robust(paper), {\sigma^2}=0.01','Robust(paper), {\sigma^2}=0.002','Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002')



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
semilogy(tx_power_arr,avg_ser_robust(1,:),'-o',tx_power_arr,avg_ser_robust(2,:),'-square',tx_power_arr,avg_ser_robust(3,:),'-*',tx_power_arr,theo_err_prb,'--');
colororder([1 0 0;0 0 1;0 0 0;1 0 1]);
xlabel("Transmit Power (dBm)")
ylabel("SER")
legend('Robust, {\sigma^2}=0.05','Robust, {\sigma^2}=0.01','Robust, {\sigma^2}=0.002','Theoretical')
