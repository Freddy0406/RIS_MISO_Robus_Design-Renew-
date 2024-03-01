clc;
clear;
close all;


% 定义参数
M = 16;      % 符号数
symbols = randi([0 M-1], 1000, 1);  % 生成随机符号序列

% 16QAM 调制
modulated_signal = qammod(symbols, M);

% 绘制调制后的信号点图
scatterplot(modulated_signal);
title('16QAM Modulated Signal');

% 16QAM 解调
demodulated_symbols = qamdemod(modulated_signal, M);


% 计算符号错误率
symbol_error_rate = sum(demodulated_symbols ~= symbols) / length(symbols);

disp(['符号错误率: ' num2str(symbol_error_rate)]);