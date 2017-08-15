%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation Interference Alignment for uncoded SC-FDMA transmission    %
%   MIMO Strictly Linear Transceiver Design for Interference Alignment    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Inputs                                %
% bits_per_symb:  how many bits per symbol for each user                  %
% offset       :  offset of each channel for each user                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; close all;
clc;

%% Variable initialization

% SNR
EbNo      = 30; % in dB
EbNo_lin  = 10.^(EbNo./10);

% User specific parameters
bits_per_symb  = [2 2];
N_user         = length(bits_per_symb); % Number of transmitter-receiver pairs

% Channel parameters
n_ch        = 1;
type        = 'A';
N_snapshot  = 20000;
Nr          = 2;
Nt          = 2;

start       = 1;
offset      = [0 200 434] + start - 1; % starting point for

% f = generate_channel(type, Nr, Nt, N_snapshot);
load([ 'LTE_channel_' type '_Anz' num2str(N_snapshot) '_cell_NR' num2str(Nr) '_NT' num2str(Nt) '.mat' ]);

% SC-FDMA specific parameters
Scfdma.l_cp      = 144;
Scfdma.M         = 300;
Scfdma.N         = 512;
Scfdma.nu_0      = 60;
Scfdma.N_scSymb  = 14;

N_x = Scfdma.N + Scfdma.l_cp;

% WL Filter Input Parameters
WL.B  = 1;                 % number of symbols per streams 
WL.br = 0;                 % number of real valued symbols
WL.Br = 2*WL.B - WL.br;
WL.nr = 0;                 % real valued components at input of each antenna
WL.Dr = 2*Nt - WL.nr;

WL.Rs = zeros(WL.Br, WL.Br, N_user);

% Get frame parameters, noise variance, BER
idx_ch    = zeros(N_user, N_snapshot);
num_bits  = zeros(1,N_user); % N bits to be generated
VarN      = zeros(N_user, length(EbNo_lin));

BER       = zeros(N_user, length(EbNo));

for i_user = 1:N_user
    % channel index vector
    i_temp = 1:N_snapshot;
    idx_ch(i_user,:) = circshift(i_temp,[0 offset(i_user)]);
    
    % frame parameters
    num_bits(i_user)     =  Scfdma.M* bits_per_symb(i_user)*Scfdma.N_scSymb;
end

for i_user = 1:N_user
    VarN(i_user,:)    = 1./(bits_per_symb(i_user).*EbNo_lin);
    if bits_per_symb(i_user) == 1           % Real valued symbols
        WL.Rs(:, :, i_user) = eye(WL.Br);
    elseif bits_per_symb(i_user) == 2       % Complex valued symbols
        WL.Rs(:, :, i_user) = 0.5*eye(WL.Br);
    end
end

tx_bits   = cell(N_user,1);
tx        = cell(N_user,1);
SymbTab   = cell(N_user,1);


for i_user = 1:N_user
    tx_bits{i_user} = randi([0 1], 1, num_bits(i_user));
    % Bit mapping
    [tx{i_user}, SymbTab{i_user}] = bitMap(tx_bits{i_user}, bits_per_symb(i_user));
    for nt = 1:Nt
        tx_bits{i_user}(nt,:)  = tx_bits{i_user};
        tx{i_user}(nt,:)       = tx{i_user};
    end
end

tx_sc  = cell(N_user,1);
count  = 1;
% check = cell(1, N_user);
tic
for i_ebNo = 1:length(EbNo)
    numBits    = zeros(Nr, N_user);
    bitError   = zeros(Nr, N_user);
    
    for i_ch = 1:n_ch
        %% Transmitter
        
        % Assign Channels for different users
        h_ch  = cell(N_user, N_user);
        H_ch  = cell(N_user, N_user);
        H_all = cell(N_user, N_user);
        for i_user = 1:N_user
            for j_user = 1:N_user
                h_ch{i_user, j_user} = h_cell{count};
                H_ch{i_user, j_user} = fft(h_cell{count}, Scfdma.N, 3);
                h_real = real(H_ch{i_user, j_user});
                h_imag = imag(H_ch{i_user, j_user});
                H_all{i_user, j_user} = [h_real -h_imag; h_imag h_real];
                count = count + 1;
            end
        end
        
        [V, G] = myPrecodingWL(H_all, H_ch, VarN(:, i_ebNo), Scfdma.N, WL);
        
        for i = 1:Scfdma.N
            for i_user = 1:N_user
                for j_user = 1:N_user 
                    check(i_user, j_user, i) = G{i_user}(:,:,i)' * H_ch{i_user, j_user}(:,:,i) * V{i_user}(:,:,i);
                end
            end
        end
        
%         error = (check(i_user, i_user, :) - 1)^2 + sum(check(i_user, i_user, :)^2) + norm(G{i_user}(:,:,i))^2*VarN(1);
        
        for i_user = 1:N_user
            tx_sc{i_user} = scfdma_tx(tx{i_user}, Scfdma, V{i_user});
        end
        %% Transmission
        
        rx_sc = cell(N_user, 1);
        for i_user = 1:N_user
            rx_sc{i_user} = scfdma_ch(tx_sc{i_user}, h_ch(i_user, :), VarN(i_user, i_ebNo), 'awgn');
        end
        %% Receiver
        
        rx       = cell(N_user,1);
        rx_bits  = cell(N_user,1);
        
        for i_user = 1:N_user
            
%             for n = 1:Scfdma.N
% %                 H_users = 0;
% %                 for j_user = 1:N_user
% %                     H_users = H_users + H_ch{i_user,j_user}(:,:,n)'*H_ch{i_user,j_user}(:,:,n);
% %                 end
%                 g(:,:,n) = (H_ch{i_user}(:,:,n)'*H_ch{i_user}(:,:,n) + VarN(i_user, i_ebNo)*eye(Nr))^-1 * H_ch{i_user}(:,:,n)';
%                 check(:,:,n) = g(:,:,n)*H_ch{i_user}(:,:,n);
%             end
            
            rx{i_user}       = scfdma_rx(rx_sc{i_user}, Scfdma, G{i_user});
            
            for inr = 1:Nr
                rx_bits{i_user}(inr,:)   = myDemapping(rx{i_user}(inr,:), bits_per_symb(i_user));
                bitError(inr, i_user)    = bitError(inr, i_user) + sum(xor(tx_bits{i_user}(inr,:), rx_bits{i_user}(inr,:)));
                numBits(inr, i_user)     = numBits(inr, i_user) + length(rx_bits{i_user}(inr,:));
            end
        end
    end
    if Nr == 1
        BER(:, i_ebNo) = bitError./numBits;
    else
        BER(:, i_ebNo) = sum(bitError)./sum(numBits);
    end
end
toc
%% Plotting BER
figure;
for i = 1:N_user
    semilogy(EbNo, BER(i,:), 'o-', 'Linewidth', 2);
    hold on;
end
title(['MIMO, ','Ch: "ITU-P', type, '", Snapshots = ', num2str(N_snapshot), ', Simulations = ', num2str(n_ch)]);
xlabel('Eb/No(dB)');
ylabel('BER')
grid on;
legend('User-1', 'User-2', 'User-3');