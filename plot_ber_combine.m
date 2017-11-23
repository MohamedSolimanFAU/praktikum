SNR = 0:5:20;
ber = zeros(3, length(SNR));

for n = 1:length(SNR)
    load(['results/SL_Nu3_Nt2_Nr2_SNR' num2str(SNR(n)) '.mat']);
    ber(:,n) = BER;

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