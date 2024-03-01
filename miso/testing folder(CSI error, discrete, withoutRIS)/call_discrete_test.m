clc;
clear;
close all;


N = 100;
bits = 3;

% Initialize theta (Start with random phases)
min = 0;
max = 2*pi;
v = exp(1i*(min+rand(N,1)*(max-min)));

partition = power(2,bits);
phase_array = 0:(2*pi/partition):2*pi;

[new_v] = discetet_phase(bits,v);
original = wrapTo2Pi(angle(v));
after = wrapTo2Pi(angle(new_v));
