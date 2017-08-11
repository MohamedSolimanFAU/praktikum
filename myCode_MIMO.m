clear variables; close all;
clc;

%%
% Channel
offset      = [0 200 434];
type        = 'ITU-PA';
N_snapshot  = 9000;
N_o         = 2;
N_i         = 2;

load([ 'LTE_channel_' type '_Anz' num2str(N_snapshot) '_cell_NR' num2str(1) '_NT' num2str(1) '.mat' ]);

%%
EbNo = 0:30; % in dB
EbNo_lin = 10.^(EbNo./10);
BER = zeros(1, length(EbNo));


% User
bits_per_symb = 2;

% SC-FDMA
l_cp = 144;
M    = 300;
N    = 512;
nu_0 = 60;

N_scSymb = 14;

N_x = N+l_cp;


num_bits = N_scSymb * M * bits_per_symb;
tx_bits  = randi([0 1], 1, num_bits);
tx_bits  = [tx_bits; tx_bits];


tx = [myMapping(tx_bits(1,:), bits_per_symb); myMapping(tx_bits(2,:), bits_per_symb)];
% tx = kron(tx, ones(N_i, 1));

Tx      = zeros(N_i, M*N_scSymb);
Tx_x    = zeros(N_i, N*N_scSymb);
tx_x    = zeros(N_i, N*N_scSymb);
tx_cp   = zeros(N_i, N*N_scSymb);

VarN    = 1./(bits_per_symb.*EbNo_lin);
N_user  = 1;
idx_ch  = 1:N_snapshot/(N_i*N_o);
n_ch    = 10;

% Generate MIMO Channels and its Hermitian
j_ch = 1;
h_mimo = cell(1, N_snapshot/(N_i*N_o));

for i_ch = 1:N_i*N_o:N_snapshot
    h_mimo{j_ch} = [h_cell{i_ch}(1) h_cell{i_ch+1}(1); h_cell{i_ch+2}(1) h_cell{i_ch+3}(1)];
    j_ch = j_ch +1;
end

tic
for i_ebNo = 1:length(EbNo)
    numBits = 0;
    bitError1 = 0;
    bitError2 = 0;
    for i_ch = 1:n_ch
        %% Transmitter
        tx_temp = tx;
        h_ch    = cell(1, N_user);
        h_ch{1} = h_mimo{i_ch};
        
        [v, g] = myPrecoding(h_ch, N_user, VarN(i_ebNo));
        V = fft(v{1}, M, 2);
        G = fft((g{1}.')', M, 2);
        
        for i_bl = 0:N_scSymb-1
            
            Tx(:, i_bl*M+(1:M)) = V .* fft(tx_temp(:, i_bl*M+(1:M)), M, N_i)./sqrt(M);
                        
            Tx_x(:, i_bl*N+nu_0+(1:M)) = Tx(:, i_bl*M+(1:M));
            
            tx_x(:, i_bl*N+(1:N)) = ifft(Tx_x(:, i_bl*N+(1:N)),N, N_i).*sqrt(N);
            
            tx_cp(:, i_bl*N_x+(1:N_x)) = [tx_x(:, i_bl*N+(N-l_cp+1:N)) tx_x(:, i_bl*N+(1:N))];
        end
        
        tx_sc = tx_cp;
        %% Transmission
        
        rx_fc = convLin(tx_sc,h_ch{1});
        noise = sqrt(0.5*VarN(i_ebNo)).*(randn(N_i,length(rx_fc))+ 1i*randn(N_i,length(rx_fc)));
        rx_sc = rx_fc+noise;
        
        %% Receiver
        
        Rx   = zeros(N_o, M*N_scSymb);
        Rx_x = zeros(N_o, N*N_scSymb);
        rx_x = zeros(N_o, N*N_scSymb);
        
        rx_temp = zeros(N_o,N_scSymb*M);
        
        rx_cp = rx_sc;
        
        
        % Equalizer
%         H  = h_ch{1};
%         Hh = H^-1;
%         F  = (Hh*H + VarN(i_ebNo)*eye(2))^-1 * Hh;
%         F = F./sqrt(M);
        
        for i_bl = 0:N_scSymb-1
            % Remove cyclic prefix
            rx_x(:, i_bl*N+(1:N)) = rx_cp(:, i_bl*N_x + (l_cp + 1:N_x));
            % Frequency domain: sc assignment/SC mapping
            Rx_x(:, i_bl*N + (1:N)) = fft(rx_x(:, i_bl*N + (1:N)), N, N_o)./sqrt(N);
            
            Rx(:, i_bl*M + (1:M)) = G .* Rx_x(:, i_bl*N + nu_0 +(1:M));
            % Equalization
            rx_temp(:, i_bl*M + (1:M)) = ifft(Rx(:, i_bl*M + (1:M)), M, N_o).* sqrt(M);
            
        end
        rx = rx_temp;
        rx_bits = [myDemapping(rx(1,:), bits_per_symb); myDemapping(rx(2,:), bits_per_symb)];
        
        bitError1 = bitError1 + sum(xor(tx_bits(1,:), rx_bits(1,:)));
        bitError2 = bitError2 + sum(xor(tx_bits(2,:), rx_bits(2,:)));
        numBits  = numBits + length(rx_bits(1,:));
    end
    BER(i_ebNo) = (bitError1 + bitError2)./(2*numBits);
end
toc


%% Plotting BER
figure;
semilogy(EbNo, BER);
title(['MIMO, ','Channel type: "', type, '", Total snapshots = ', num2str(N_snapshot)]);
legend(['Simulated realisations = ', num2str(n_ch)]);
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;

% figure;
% plot(rx(1,:), 'o'); hold on; plot(tx(1,:), 'o');
% 
% figure;
% plot(rx(2,:), 'o'); hold on; plot(tx(2,:), 'o');
%% TEST
% Equalization
% W   = cellfun(@squeeze,num2cell(h_ch{i_ch},3),'UniformOutput' , false);
% Wh  = cellfun(@ctranspose, W, 'UniformOutput' , false);
% Wm  = cellfun(@(x,y) x*y, Wh,W, 'UniformOutput', false);
% Wa  = cellfun(@(x,y) x+y, Wm, num2cell(VarN(i_ebNo)*eye(2)), 'UniformOutput', false);
% Wi  = cellfun(@inv, Wa, 'UniformOutput', false);
% WF  = cellfun(@(x,y) x*y, Wi,Wh, 'UniformOutput', false);
% G   = cellfun(@fft, WF, num2cell(300*ones(2)),  'UniformOutput', false);
% G   = cellfun(@(x,y) x/y, G, num2cell(sqrt(M)*ones(2)), 'UniformOutput', false);