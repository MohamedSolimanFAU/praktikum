%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation Interference Alignment for uncoded SC-FDMA transmission    %
%   SISO Strictly Linear Transceiver Design for Interference Alignment    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Inputs                                %
% bits_per_symb:  how many bits per symbol for each user                  %
% offset       :  offset of each channel for each user                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all;
clc;

%% System Configuration

% SNR
EbNo      = 0:20; % in dB
EbNo_lin  = 10.^(EbNo./10);

% User specific parameters
bits_per_symb  = [2 2 2];
N_user         = length(bits_per_symb); % Number of transmitter-receiver pairs

% SC-FDMA specific parameters
Scfdma.l_cp      = 144;
Scfdma.M         = 300;
Scfdma.N         = 512;
Scfdma.nu_0      = 60;
Scfdma.N_scSymb  = 14;

N_x = Scfdma.N + Scfdma.l_cp;

% Channel parameters
n_ch        = 2;
type        = 'ITU-PA';
N_snapshot  = 9000;
Nr          = 1;
Nt          = 1;

start       = 1;
offset      = [0 200 434] + start - 1; % starting point for


load([ 'LTE_channel_' type '_Anz' num2str(N_snapshot) '_cell_NR' num2str(Nr) '_NT' num2str(Nt) '.mat' ]);

% Get frame parameters, noise variance, BER
idx_ch    = zeros(N_user, N_snapshot);
num_bits  = zeros(1,N_user);
VarN      = zeros(N_user, length(EbNo_lin));
BER       = zeros(N_user, length(EbNo));

for i_user = 1:N_user
    % channel index vector
    i_temp = 1:N_snapshot;
    idx_ch(i_user,:) = circshift(i_temp,[0 offset(i_user)]);
    
    % frame parameters
    num_bits(i_user)     =  Scfdma.M* bits_per_symb(i_user)* Scfdma.N_scSymb;
end

for i_user = 1:N_user
    VarN(i_user,:)    = 1./(bits_per_symb(i_user).*EbNo_lin);
end


tx_bits   = cell(1,N_user);
tx        = cell(1,N_user);
SymbTab   = cell(1,N_user);

for i_user = 1:N_user
    tx_bits{i_user} = randi([0 1], 1, num_bits(i_user));
    % Bit mapping
    [tx{i_user}, SymbTab{i_user}] = bitMap(tx_bits{i_user}, bits_per_symb(i_user));
end

tx_sc   = cell(1,N_user);
M = Scfdma.M;
N = Scfdma.N;

tic
for i_ebNo = 1:length(EbNo)
    numBits   = zeros(1,N_user);
    bitError  = zeros(1,N_user);
    
    for i_ch = 1:n_ch
        %% Transmitter
                
        h =zeros(N_user, N_user);
        for i_user = 1:N_user
            temp = 1:N_user;
            idx = circshift(temp, [0 i_user-1]);
            h(i_user,:) = squeeze(h_cell{idx_ch(i_user,i_ch)}(idx))';
        end
        h_ch = num2cell(h);
        
        [v, g] = myPrecoding(h_ch, N_user, VarN(:, i_ebNo), N);
        
        for i_user = 1:N_user
            V = fft(v{i_user}, M);
            tx_sc{i_user} = scfdma_tx(tx{i_user}, Scfdma, V);
        end
        %% Transmission
        rx_sc = cell(1,N_user);
        for i_user = 1:N_user
            % rx_sc{i_user} = scfdma_ch(tx_sc{i_user}, h(i_user, :), VarN(i_user, i_ebNo), 'noiseless');
            rx_sc{i_user} = scfdma_ch(tx_sc{i_user}, h_ch{i_user, i_user}, VarN(i_user, i_ebNo), 'awgn');
        end
        %% Receiver
        
        rx       = cell(N_user,1);
        rx_bits  = cell(N_user,1);
        
        for i_user = 1:N_user
            G = fft(g{i_user}', M);
            
            rx{i_user}       = scfdma_rx(rx_sc{i_user}, Scfdma, G);
            
            rx_bits{i_user}  = myDemapping(rx{i_user}, bits_per_symb(i_user));
            
            bitError(i_user) = bitError(i_user) + sum(xor(tx_bits{i_user}, rx_bits{i_user}));
            numBits(i_user)  = numBits(i_user) + length(rx_bits{i_user});
        end
    end
    BER(:, i_ebNo) = bitError./numBits;
end
toc


%% Plotting BER
figure;
for i = 1:N_user
    semilogy(EbNo, BER(i,:), 'o-', 'Linewidth', 2);
    hold on;
end
title(['SISO, ','Ch: "', type, '", Snapshots = ', num2str(N_snapshot), ', Simulations = ', num2str(n_ch)]);
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
legend('User-1', 'User-2', 'User-3');

% semilogy(EbNo, BER(1,:), 'bo-', 'Linewidth', 2);
% title(['SISO, ','Ch: "', type, '", Snapshots = ', num2str(N_snapshot), ', Simulations = ', num2str(n_ch)]);
% xlabel('Eb/No(dB)');
% ylabel('BER')
% grid on;
% hold on;
%
% semilogy(EbNo, BER(2,:), 'ro-', 'Linewidth', 2);
% semilogy(EbNo, BER(3,:), 'go-', 'Linewidth', 2);