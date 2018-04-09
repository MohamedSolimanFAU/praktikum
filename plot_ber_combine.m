clc; clear all; close all;

Nu = 3;

SNR = 0:5:25;
ber  = zeros(Nu, length(SNR));
bler = zeros(Nu, length(SNR));

for n = 1:length(SNR)
    load(['results/cluster/SL_Nu3_N22_BPSK_PA_SNR_' num2str(SNR(n)) '.mat']);
    ber(:,n)  = BER;
    bler(:,n) = BLER;
end

figure;
for i = 1:N_user
    semilogy(SNR, ber(i,:), 'o-', 'Linewidth', 2);
    hold on;
end
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
legend('User-1', 'User-2', 'User-3');
title(['ITU-PA, ', 'BPSK ', 'Nusers = 3']);


ber  = zeros(Nu, length(SNR));
bler = zeros(Nu, length(SNR));
for n = 1:length(SNR)
    load(['results/cluster/WL_Nu4_N22_BPSK_PA_SNR_' num2str(SNR(n)) '.mat']);
    ber(:,n)  = BER;
    bler(:,n) = BLER;
end

for i = 1:N_user
    semilogy(SNR, ber(i,:), 'o--', 'Linewidth', 2);
    hold on;
end
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
legend('SL-User1', 'SL-User2', 'SL-User3', 'WL-User1', 'WL-User2', 'WL-User3');
title(['ITU-PA, ', 'BER for BPSK', ', Nusers = 3']);

% figure;
% for i = 1:N_user
%     semilogy(SNR, bler(i,:)./800, 'o-', 'Linewidth', 2);
%     hold on;
% end
% xlabel('Eb/No(dB)');
% ylabel('BER')
% grid on;
% legend('User-1', 'User-2', 'User-3', 'User-4', 'User-5');
% title(['ITU-PA, ', 'Widely linear filter ', 'Nusers = 3']);