clc;
clear
close all;

%% Initialization parameter
    N = 40;
    P0 = 1;
    M = 4;                          %BS antenna amount
    variance = 0.01;
    P0_antenna = P0/M;
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
    power_array = zeros(1,10000);

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


    
    % Emulate estimation channel

    K = 10;                                     %Rician channel factor
    G = zeros(N,M);                             %Ideal channel => Rician fading
    SqrtKoverKplus1 = sqrt(K / (K + 1));
    Sqrt1overKplus1 = sqrt(1 / (K + 1));
    Sqrt2 = sqrt(2);
    SqrtG_PL = sqrt(G_PL);
    Sqrthr_PL = sqrt(hr_PL);
    Sqrthd_PL = sqrt(hd_PL);
    delta_coef = sqrt(variance)/Sqrt2;

    for a = 1:N
        for b = 1:M
            los =  SqrtKoverKplus1*randn(1) + 1i* SqrtKoverKplus1*randn(1);
            multipath = Sqrt1overKplus1*randn(1) + 1i*Sqrt1overKplus1*randn(1);
            G(a,b) = (los+multipath)/Sqrt2;                %Ideal channel => Rician fading
        end
    end
    % Let norm(G)^2 = G_PL
    
    for i = 1:M
        G(:,i) = SqrtG_PL*G(:,i)/norm(G(:,i));
    end
    

    hr = zeros(N,1);
    for num = 1:N
        los = SqrtKoverKplus1*randn(1) + 1i*SqrtKoverKplus1*randn(1);
        multipath = Sqrt1overKplus1*randn(1) + 1i*Sqrt1overKplus1*randn(1);
        hr(num) = (los+multipath)/Sqrt2;                %Ideal channel => Rician fading
    end
    % Let norm(hr)^2 = hr_PL
    hr = Sqrthr_PL*hr/norm(hr);

    hd = ((randn(M,1)+1i*randn(M,1)))/Sqrt2;      %Ideal channel => Rayleigh-flat fading (Gaussian)                         
    % Let norm(hd)^2 = hd_PL
    hd = Sqrthd_PL*hd/norm(hd);

    G_hat = G;         %G = G_hat+delta_G ; G_hat=G-delta_G        (AP-IRS)
    hd_hat = hd;       %hd = hd_hat+delta_hd ; hd_hat=hd-delta_hd  (IRS-User)
    hr_hat = hr;       %hr = hr_hat+delta_hr ; hr_hat=hr-delta_hr  (AP-User)

    % Generate delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
    % Let norm(delta_G)^2 = G_PL ; norm(delta_hd)^2 = hd_PL ; norm(delta_hr)^2 = hr_PL
    delta_G = (randn(N,M)+1i*randn(N,M))*delta_coef;
    delta_G = SqrtG_PL*delta_G/norm(delta_G);
    delta_hd = (randn(M,1)+1i*randn(M,1))*delta_coef;
    delta_hd = Sqrthd_PL*delta_hd*norm(delta_hd); 
    delta_hr = (randn(N,1)+1i*randn(N,1))*delta_coef;
    delta_hr = Sqrthr_PL*delta_hr/norm(delta_hr);







    %% SER calculation
%     count_of_SYM = 0;
%     numoferr=0;
%     numofbiterr=0;
% 
%     MOrder=16;
%     N0 = noise_var;
%     Eav= P0;		 	  	        % energy per symbol
%     symbol_d= sqrt(P0/10);	        % min. distance between symbols
%     qam_sig = zeros(1,2);
%     de_sig =  zeros(1,2);
%     r = zeros(1,2);
% 
%     % Mapping to the signal constellation follows.
%     mapping=[-3*symbol_d 3*symbol_d;
%               -symbol_d  3*symbol_d;
%                symbol_d  3*symbol_d;
%              3*symbol_d  3*symbol_d;
%               -3*symbol_d  symbol_d;
%                  -symbol_d  symbol_d;
%                  symbol_d  symbol_d;
%               3*symbol_d  symbol_d;
%                  -3*symbol_d  -symbol_d; 
%                -symbol_d  -symbol_d; 
%             symbol_d  -symbol_d;
%           3*symbol_d  -symbol_d;
%          -3*symbol_d  -3*symbol_d;
%            -symbol_d  -3*symbol_d;
%             symbol_d  -3*symbol_d;
%           3*symbol_d  -3*symbol_d];
% 
%         temp=rand;		        	  	% a uniform R.V. between 0 and 1
%         dsource=1+floor(MOrder*temp);	  	% a number between 1 and 16, uniform  
%         qam_sig=mapping(dsource,:);
%         [dsource1,dsource2,dsource3,dsource4] = determin_bit(qam_sig,symbol_d);
% 
% 
%         % received signal
%         SqrtNoise_var = sqrt(noise_var);
%         n  = [SqrtNoise_var*randn SqrtNoise_var*randn];
%         noise = (n(1)+1i*n(2))/Sqrt2;
%         s = (qam_sig(1,1)+1i*qam_sig(1,2))/Sqrt2;
%         y = (hr'*btheta*G+hd')*w*s+noise;



