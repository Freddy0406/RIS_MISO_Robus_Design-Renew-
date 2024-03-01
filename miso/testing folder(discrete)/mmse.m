function [mse,avg_SER,power_array]=mmse(N,variance,P0,mode,bit_of_phase)
%% Initialization parameter
    %N  Reflecting element amount per RIS
    %P0 Maximum transmit power of BS
    %variance var_g=var_r=var_d=var
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

    % Generate delta_G, delta_hd, delta_hr (AP-IRS, IRS-User, AP-User link)
    % Let norm(delta_G)^2 = G_PL ; norm(delta_hd)^2 = hd_PL ; norm(delta_hr)^2 = hr_PL
    delta_G = (randn(N,M)+1i*randn(N,M))*delta_coef;
    delta_G = SqrtG_PL*delta_G/norm(delta_G);
    delta_hd = (randn(M,1)+1i*randn(M,1))*delta_coef;
    delta_hd = Sqrthd_PL*delta_hd/norm(delta_hd);
    delta_hr = (randn(N,1)+1i*randn(N,1))*delta_coef;
    delta_hr = Sqrthr_PL*delta_hr/norm(delta_hr);

    for a = 1:N
        for b = 1:M
            los =  SqrtKoverKplus1*randn(1) + 1i* SqrtKoverKplus1*randn(1);
            multipath = Sqrt1overKplus1*randn(1) + 1i*Sqrt1overKplus1*randn(1);
            G(a,b) = (los+multipath)/Sqrt2;                %Ideal channel => Rician fading
        end
    end
    % Let norm(G)^2 = G_PL
    G = SqrtG_PL*G/norm(G);

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

    G_hat = G-delta_G;          %G = G_hat+delta_G ; G_hat=G-delta_G        (AP-IRS)
    hd_hat = hd-delta_hd;       %hd = hd_hat+delta_hd ; hd_hat=hd-delta_hd  (IRS-User)
    hr_hat = hr-delta_hr;       %hr = hr_hat+delta_hr ; hr_hat=hr-delta_hr  (AP-User)
    mse_before = 1.0;


    
    %% Start optimization
    variance_G = variance*G_PL;
    variance_hd = variance*hd_PL;
    variance_hr = variance*hr_PL;
    t=1;                        %Set times
%% mode 1:The proposed robust design
    if (mode==1)           
        while true
            fprintf('The %d generation...\n',t);
            %% Optimization of c
            fprintf('\tOptimization of c...\n');        
            A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                hd_hat*hd_hat'+variance_G*power(norm(hr_hat),2)*eye(M)+(variance_hr*G_hat')*G_hat+...
                (N*variance_hr*variance_G+variance_hd)*eye(M);

            alpha = (G_hat'*btheta'*hr_hat+hd_hat);
        
            c = (w'*alpha)/(w'*A*w+noise_var);          %Wiener Filter
        
            %% Optimization of w
            fprintf('\tOptimization of w...\n');
            lambda = 0.0;
            w =  pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));  %Assume lambda(Lagrange mult.) = 0

            if(power(norm(w),2)>P0)
                [lambda] = search(c,alpha,A,P0,M);
                w = pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));
            else
            end
            power_array(t) = power(norm(w),2);
            
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
            if(abs(mse_before-mse)>epsilon)
                mse_before = mse;
            elseif(abs(mse_before-mse)<=epsilon)
                [final_v] = discetet_phase(bit_of_phase,v);
                v = final_v;
                btheta = diag(v);
                A = (G_hat'*btheta'*hr_hat)*hr_hat'*btheta*G_hat+...          %According to the formula              
                    hd_hat*hr_hat'*btheta*G_hat+G_hat'*btheta'*hr_hat*hd_hat'+...
                    hd_hat*hd_hat'+variance_G*power(norm(hr_hat),2)*eye(M)+(variance_hr*G_hat')*G_hat+...
                    (N*variance_hr*variance_G+variance_hd)*eye(M);

                alpha = (G_hat'*btheta'*hr_hat+hd_hat);
                mse = power(abs(c),2)*(w'*A*w+noise_var)-w'*alpha*conj(c)-c*alpha'*w+1;
                break;
            end             
            t = t+1;
        end
        %% SER
%         symbol_error_rate = zeros(1,100);
%         M = 16;                                 % M-QAM
%         modulated_signal = zeros(1000,1);
%         y = zeros(size(modulated_signal));
%         y_de = zeros(size(modulated_signal));
%         s_hat = zeros(size(modulated_signal));
%         for j = 1:100
%             
%             symbols = randi([0 M-1], 1000, 1);     % Random integer between 0~M
%             modulated_signal = qammod(symbols, M);  % M-QAM modulation
%             
%             for i = 1:length(modulated_signal)
%                 s = modulated_signal(i);
%                 noise = sqrt(noise_var)*randn(1,1,"like",1i);   %Generate random 1X1 noise
%                 y(i) = (hr_hat'*btheta*G_hat+hd_hat')*w*s+noise;
%                 s_hat(i) = c*y(i);
%                 y_de(i) = qamdemod(s_hat(i), M);
%             end
%             symbol_error_rate(j) = sum(y_de ~= symbols) / length(symbols);
%         end
%         avg_SER = mean(symbol_error_rate);
            avg_SER = 0.0;
    end
end
