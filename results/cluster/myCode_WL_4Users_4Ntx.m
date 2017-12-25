%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Simulation Interference Alignment for uncoded SC-FDMA transmission    %
%   MIMO Strictly Linear Transceiver Design for Interference Alignment    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Inputs                                %
% bits_per_symb:  how many bits per symbol for each user                  %
% offset       :  offset of each channel for each user                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Variable initialization

EbNo_lin  = 10.^(EbNo./10);

% User specific parameters
bits_per_symb  = [2 2 2 2];
N_user         = length(bits_per_symb); % Number of transmitter-receiver pairs

% Channel parameters
n_ch        = 1000;
type        = 'ITU-PA';
N_snapshot  = 20000;
Nr          = 2;
Nt          = 4;

start       = 1;
offset      = [0 200 434 656] + start - 1; % starting point for

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

numSubframes    = zeros(N_user, length(EbNo));
subframeError   = zeros(N_user, length(EbNo));
new_bitError    = zeros(N_user,1);
BLER            = zeros(N_user, length(EbNo));

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

tx_bits    = cell(N_user, 1);
tx_symbols = cell(N_user, 1);
tx         = cell(N_user, 1);
SymbTab    = cell(N_user, 1);


for i_user = 1:N_user
    tx_bits{i_user} = randi([0 1], 1, num_bits(i_user));
    % Bit mapping
    [tx{i_user}, SymbTab{i_user}] = bitMap(tx_bits{i_user}, bits_per_symb(i_user));
end

tx_sc  = cell(N_user, 1);
check  = cell(N_user, 1);

tic
for i_ebNo = 1:length(EbNo)
    count  = 1;
    
    numBits    = zeros(1, N_user);
    bitError   = zeros(1, N_user);
    
    numSubframes(:, i_ebNo) = numSubframes(:, i_ebNo) + ones(i_user,1);
    
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
        
        % apply precoding
        [V, G] = myPrecodingWL(H_all, H_ch, VarN(:, i_ebNo), Scfdma.N, WL);
        
        %H_ov
        for i = 1:Scfdma.N
            for i_user = 1:N_user
                for j_user = 1:N_user 
                    check{i_user, j_user}(:,:,i) = G{i_user}(:,:,i)' * H_all{i_user, j_user}(:,:,i) * V{i_user}(:,:,i);
                end
            end
        end

        % transmitter
        for i_user = 1:N_user
            tx_sc{i_user} = scfdma_txWL(tx{i_user}, Scfdma, V{i_user});
        end
        %% Transmission
        
        rx_sc = cell(N_user, 1);
        for i_user = 1:N_user
            rx_sc{i_user} = scfdma_ch(tx_sc, h_ch(i_user, :), VarN(i_user, i_ebNo), 'awgn');
        end
        %% Receiver
        
        rx       = cell(N_user,1);
        rx_bits  = cell(N_user,1);
        
        for i_user = 1:N_user
            rx{i_user}       = scfdma_rxWL(rx_sc{i_user}, Scfdma, G{i_user});
            rx_bits{i_user}  = myDemapping(rx{i_user}, bits_per_symb(i_user));
            
            bitError(1, i_user)    = bitError(1, i_user) + sum(xor(tx_bits{i_user}, rx_bits{i_user}));
            numBits(1, i_user)     = numBits(1, i_user) + length(rx_bits{i_user});
            
            new_bitError(i_user,1) = sum(tx_bits{i_user} ~= rx_bits{i_user});
        end
        subframeError_new        = new_bitError > 0;
        subframeError(:, i_ebNo) = subframeError(:, i_ebNo) + subframeError_new;

    end
    BER(:, i_ebNo)  = bitError./numBits;
    BLER(:,i_ebNo)  = subframeError(:,i_ebNo)./numSubframes(:,i_ebNo);

    mat_file_name = strcat('results/WL_Nu4_Nt4_Nr2_PA_SNR_', num2str(EbNo(i_ebNo)), '.mat');
    save(mat_file_name , 'BER', 'BLER', 'bitError', 'N_user', 'n_ch', 'Nr', 'Nt', 'subframeError', 'numSubframes', 'i_ch');
end
toc
