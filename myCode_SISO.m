clear variables; close all;
clc;

%%
% Channel
offset      = [0 200 434];
type        = 'ITU-PA';
N_snapshot  = 9000;
N_o         = 1;
N_i         = 1;

load([ 'LTE_channel_' type '_Anz' num2str(N_snapshot) '_cell_NR' num2str(N_o) '_NT' num2str(N_i) '.mat' ]);

% clear N_o; clear N_i;

%%
EbNo = 0:30; % in dB
EbNo_lin = 10.^(EbNo./10);
BER = zeros(1, length(EbNo));


% User
bits_per_symb = 2;% [2 2 2];

% SC-FDMA
l_cp = 144;
M    = 300;
N    = 512;
nu_0 = 60;

N_scSymb = 14;

N_x = N + l_cp;


num_bits = N_scSymb * M * bits_per_symb;
tx_bits  = randi([0 1], 1, num_bits);

tx = myMapping(tx_bits, bits_per_symb);


Tx      = zeros(1, M*N_scSymb);
Tx_x    = zeros(1, N*N_scSymb);
tx_x    = zeros(1, N*N_scSymb);
tx_cp   = zeros(1, N*N_scSymb);

VarN    = 1./(bits_per_symb.*EbNo_lin);
N_user  = 1;
idx_ch  = 1:N_snapshot;
n_ch    = 1000;

tic
for i_ebNo = 1:length(EbNo)
    numBits = 0;
    bitError = 0;
    for i_ch = 1:n_ch
        %% Transmitter
        tx_temp = tx;
        h_ch = cell(1, N_user);
        h_ch{1} = h_cell{i_ch};
        
        [v, g] = myPrecoding(h_ch, N_user, VarN(i_ebNo));
        
        for i_bl = 0:N_scSymb-1
            Tx(i_bl*M+(1:M)) = fft(v{1} .* tx_temp(i_bl*M+(1:M)),M)./sqrt(M);% 
            
            Tx_x(i_bl*N+nu_0+(1:M)) = Tx(i_bl*M+(1:M));
            
            tx_x(i_bl*N+(1:N)) = ifft(Tx_x(i_bl*N+(1:N)),N).*sqrt(N);
            
            tx_cp(i_bl*N_x+(1:N_x)) = [tx_x(i_bl*N+(N-l_cp+1:N)) tx_x(i_bl*N+(1:N))];
        end
        
        tx_sc = tx_cp;
        %% Transmission
        
        rx_fc = convLin(tx_sc,h_ch{1});
        noise = sqrt(0.5*VarN(i_ebNo)).*(randn(1,length(rx_fc))+ 1i*randn(1,length(rx_fc)));
        rx_sc = rx_fc+noise;
        
        %% Receiver
        
        Rx = zeros(1,M*N_scSymb);
        Rx_x = zeros(1,N*N_scSymb);
        rx_x = zeros(1,N*N_scSymb);
        
        rx_temp = zeros(1,N_scSymb*M);
        
        rx_cp = rx_sc;
        % MMSE Equalization filter
%         h      = squeeze(h_ch);
%         H(:,1) = fft(h, M)./sqrt(M);
%         F      = (H'*H + VarN(i_ebNo))^-1 * H';
%         F     = 1./H;
        G = fft(g{1}', M);
        
        for i_bl = 0:N_scSymb-1
            rx_x(i_bl*N+(1:N)) = rx_cp(i_bl*N_x+(l_cp+1:N_x));
            % Frequency domain: sc assignment/SC mapping
            Rx_x(i_bl*N+(1:N)) = fft(rx_x(i_bl*N+(1:N)), N)./sqrt(N);
            
            Rx(i_bl*M+(1:M)) = G .* Rx_x(i_bl*N + nu_0+(1:M));
            
            % Equalization
            rx_temp(i_bl*M+(1:M)) = ifft(Rx(i_bl*M+(1:M)), M).*sqrt(M);
            
        end
        rx = rx_temp;
        [rx_bits] = myDemapping(rx, bits_per_symb);
        
        bitError = bitError + sum(xor(tx_bits, rx_bits));
        numBits  = numBits + length(rx_bits);
    end
    BER(i_ebNo) = bitError./numBits;
end
toc


%% Plotting BER
figure;
semilogy(EbNo, BER);
title(['SISO, MMSE, ','Channel type: "', type, '", Total snapshots = ', num2str(N_snapshot)]);
legend(['Simulated realisations = ', num2str(n_ch)]);
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
% figure;
% plot(rx, 'o'); hold on; plot(tx, 'o');
% grid on

% saveas(gcf, 'filename.jpg');