clc;
clear;
close all;


%% Initialization parameter
N = 40;
P0 = round(power(10,0),5);
variance = 0.08;
M = 4;                          %BS antenna amount
L0 = power(10,(-30/10));        %Path loss dB -30dB         
AP_loc = [0 0];                 %AP location
IRS_loc = [100,0];              %IRS location
UE_loc = [100 20];              %User location

AP_UE_dis = sqrt(sum((UE_loc-AP_loc).^2));      %AP-UE distance  
AP_IRS_dis = sqrt(sum((AP_loc-IRS_loc).^2));    %AP-IRS distance
UE_IRS_dis = sqrt(sum((UE_loc-IRS_loc).^2));    %UE-IRS distance
alpha_LoS = -2;                                 %Path Loss Exponent of LoS(AP->IRS->UE)
alpha_nLoS = -3;                                %Path Loss Exponent of non-LoS(AP->UE)
G_PL =  L0*power((AP_IRS_dis),alpha_LoS);       %Path Loss of LoS(AP->IRS->UE)
hr_PL = L0*power((UE_IRS_dis),alpha_LoS);       %Path Loss of LoS(AP->IRS->UE)
hd_PL = L0*power((AP_UE_dis),alpha_nLoS);       %Path Loss of non-LoS(AP->UE)
epsilon = power(10,-4);                         %Limit of optimization iteration


%Noise
noise_var = power(10,-11);                      %var_n = 110dBm = 10^-11 mW

% Initialize w 
w = zeros(M,1);
for i = 1:M
    w(i) = sqrt(P0)/sqrt(M);
end

% Initialize theta (Start with random phases)
min_phase = 0;
max_phase = 2*pi;
v = exp(1i*(min_phase+rand(N,1)*(max_phase-min_phase)));
btheta = diag(v);

% Generate delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
delta_G = (sqrt(variance)*randn(N,M,"like",1i))/sqrt(2)/sqrt(N)*sqrt(G_PL);
delta_hd = (sqrt(variance)*randn(M,1,"like",1i))/sqrt(2)/sqrt(M)*sqrt(hd_PL);
delta_hr = (sqrt(variance)*randn(N,1,"like",1i))/sqrt(2)/sqrt(N)*sqrt(hr_PL);

% Emulate estimation channel

K = 10;                                     %Rician channel factor
G_hat = zeros(N,M);                             %Ideal channel => Rician fading
for a = 1:N
    for b = 1:M
        los = sqrt(K / (K + 1))*randn(1) + 1i*sqrt(K / (K + 1))*randn(1);
        multipath = sqrt(1 / (K + 1)) *randn(1) + 1i*sqrt(1 / (K + 1))*randn(1);
        G_hat(a,b) = (los+multipath)/sqrt(2)/sqrt(N)*sqrt(G_PL);                %Ideal channel => Rician fading           
    end
end

hr_hat = zeros(N,1);
for num = 1:N
    los = sqrt(K / (K + 1))*randn(1) + 1i*sqrt(K / (K + 1))*randn(1);
    multipath = sqrt(1 / (K + 1)) *randn(1) + 1i*sqrt(1 / (K + 1))*randn(1);
    hr_hat(num) = (los+multipath)/sqrt(2)/sqrt(N)*sqrt(hr_PL);                %Ideal channel => Rician fading
end


hd_hat = ((randn(M,1)+1i*randn(M,1))/sqrt(2))/sqrt(M)*sqrt(hd_PL);      %Ideal channel => Rayleigh-flat fading (Gaussian)                         

G = G_hat+delta_G;
hd = hd_hat+delta_hd;
hr = hr_hat+delta_hr;
%G = G_hat+delta_G ; G_hat=G-delta_G        (AP-IRS)
%hd = hd_hat+delta_hd ; hd_hat=hd-delta_hd  (IRS-User)
%hr = hr_hat+delta_hr ; hr_hat=hr-delta_hr  (AP-User)
mse_before = 0;

t = 1;
while true
    fprintf('The %d generation...\n',t);
    %% Optimization of c
    fprintf('\tOptimization of c...\n');        
    A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
        hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
        hd_hat*hd_hat'+variance*power(norm(hr_hat),2)*eye(M)+(variance*G_hat')*G_hat+...
        (N*variance*variance+variance)*eye(M);

    alpha = (G_hat'*btheta'*hr_hat+hd_hat);

    c = (w'*alpha)/(w'*A*w+noise_var);          %Wiener Filter

    %% Optimization of w
    fprintf('\tOptimization of w...\n');
    lambda = 0.0;
    w =  pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0

    if(power(norm(w),2)>P0)
        [lambda] = search(c,alpha,A,P0,M);
        w = pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));
    end
    
    %% Optimization of theta
    fprintf('\tOptimization of btheta...\n');
    Phi = diag(hr_hat')*G_hat*w*c;
    d = hd_hat'*w*c;
    Q = Phi*Phi';
    q = Phi*(1-conj(d));
    [new_v] = opt_v(Q,q,v);
    v = new_v;
    btheta = diag(v);

    %% Calculate mse
    fprintf('\tCalculate mmse...\n\n');
    mse = power(abs(c),2)*(w'*A*w+noise_var)-w'*alpha*conj(c)-c*alpha'*w+1;
    if((mse_before-mse)>epsilon)
        mse_before = mse;
    elseif((mse_before-mse)<=epsilon)
        mse = power(abs(c),2)*(w'*A*w+noise_var)-w'*alpha*conj(c)-c*alpha'*w+1;
        break;
    end             
    t = t+1;
end



    