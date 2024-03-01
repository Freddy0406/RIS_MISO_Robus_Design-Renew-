clc;
clear;
close all;


P0 = round(power(10,0),5);
noise_var = power(10,-11);
Num_of_symbol=100;
MOrder=16;
N0 = noise_var;
Eav= P0;		 	  	        % energy per symbol
symbol_d= sqrt(P0/10);	        % min. distance between symbols

dsource1 = zeros(1,Num_of_symbol);          %第1個bit
dsource2 = zeros(1,Num_of_symbol);          %第2個bit
dsource3 = zeros(1,Num_of_symbol);          %第3個bit
dsource4 = zeros(1,Num_of_symbol);          %第4個bit

dsource = zeros(1,Num_of_symbol);
for i=1:Num_of_symbol
  temp=rand;		        	  	% a uniform R.V. between 0 and 1
  dsource(i)=1+floor(MOrder*temp);	  	% a number between 1 and 16, uniform 
end
% Mapping to the signal constellation follows.
mapping=[-3*symbol_d 3*symbol_d;
          -symbol_d  3*symbol_d;
           symbol_d  3*symbol_d;
         3*symbol_d  3*symbol_d;
          -3*symbol_d  symbol_d;
             -symbol_d  symbol_d;
             symbol_d  symbol_d;
          3*symbol_d  symbol_d;
	         -3*symbol_d  -symbol_d; 
           -symbol_d  -symbol_d; 
        symbol_d  -symbol_d;
      3*symbol_d  -symbol_d;
     -3*symbol_d  -3*symbol_d;
       -symbol_d  -3*symbol_d;
        symbol_d  -3*symbol_d;
      3*symbol_d  -3*symbol_d];

for i=1:Num_of_symbol
    qam_sig(i,:)=mapping(dsource(i),:);
    [dsource1(i),dsource2(i),dsource3(i),dsource4(i)] = determin_bit(qam_sig(i,:),symbol_d);
end

% received signal
SqrtNoise_var = sqrt(noise_var);
for i=1:Num_of_symbol
  n  = [SqrtNoise_var*randn SqrtNoise_var*randn];
  noise = n(1)+1i*n(2);
  s = qam_sig(i,1)+1i*qam_sig(i,2);
  y = (hr'*btheta*G+hd')*w*s+noise;
  s_hat = c*y;
  r(i,:) = [real(s_hat) imag(s_hat)];
end

numoferr=0;
numofbiterr=0;
distance = zeros(1,16);
for i=1:Num_of_symbol
    for j = 1:16
        target = r(i,:)-mapping(j,:);
        distance(j) = target(1)^2+target(2)^2;
    end
    a= min(distance);
    decis = find(distance==a);
    de_sig(i,:)=mapping(decis,:);
    [de_1,de_2,de_3,de_4] = determin_bit(de_sig(i,:),symbol_d);

    if (de_1~=dsource1(i))
        numofbiterr_a = numofbiterr_a+1;
    elseif (de_2~=dsource2(i))
        numofbiterr_a = numofbiterr_a+1;
    elseif (de_3~=dsource3(i))
        numofbiterr_a = numofbiterr_a+1;
    elseif (de_4~=dsource4(i))
        numofbiterr_a = numofbiterr_a+1;    
    end

    if (decis~=dsource(i))
        numoferr=numoferr+1;
    end
end
avg_SER=numoferr/(Num_of_symbol);
pb = numofbiterr/(4*Num_of_symbol);
