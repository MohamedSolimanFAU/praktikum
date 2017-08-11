clear variables; close all;
clc;

%%

EbNo = 0:20; % in dB

BER_9000  = struct2array(load('BER_SISO/BER_SISO_9000.mat'));
BER_1000  = struct2array(load('BER_SISO/BER_SISO_1000.mat'));
BER_90000 = struct2array(load('BER_SISO/BER_SISO_90000.mat'));
BER_25000 = struct2array(load('BER_SISO/BER_SISO_25000.mat'));


figure;
semilogy(EbNo, BER_9000, 'b', 'Linewidth', 2);
title('SISO, N_{users} = 1, MMSE Equalization + Precoding');
xlabel('Eb/No(dB)');
ylabel('BER');
grid on;

hold on;
semilogy(EbNo, BER_90000, 'r', 'Linewidth', 2);
semilogy(EbNo, BER_25000, 'k', 'Linewidth', 2);
semilogy(EbNo, BER_1000, 'g', 'Linewidth', 2);

legend('PA 9000', 'PA 90000', 'PA 25000', 'PB 1000');

