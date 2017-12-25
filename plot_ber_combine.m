SNR = 5:5:25;
ber  = zeros(3, length(SNR));
bler = zeros(3, length(SNR));

for n = 1:length(SNR)
    load(['results/cluster/WL_Nu3_Nt2_Nr2_PA_SNR_' num2str(SNR(n)) '.mat']);
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

figure;
for i = 1:N_user
    semilogy(SNR, bler(i,:), 'o-', 'Linewidth', 2);
    hold on;
end
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
legend('User-1', 'User-2', 'User-3');