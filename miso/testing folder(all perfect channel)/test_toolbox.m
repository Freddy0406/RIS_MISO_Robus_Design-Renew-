clc;
clear;
close all;

N = 2000;
K = 10;
M = 1000;

K = 10;                                     %Rician channel factor
G = zeros(N,M);                             %Ideal channel => Rician fading
for a = 1:N
    for b = 1:M
        los = sqrt(K / (K + 1))*randn(1) + 1i*sqrt(K / (K + 1))*randn(1);
        multipath = sqrt(1 / (K + 1)) *randn(1) + 1i*sqrt(1 / (K + 1))*randn(1);
        G(a,b) = (los+multipath)/sqrt(2);                %Ideal channel => Rician fading           
    end
end
% norm(G)
% G = G/norm(G);
% norm(G)
% G = 0.05*G;
% norm(G)
% 
% s = rand(1,2000)+1i*rand(1,2000);
% s = s/norm(s);



