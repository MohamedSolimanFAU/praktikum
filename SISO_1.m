% Simulation Interference Alignment for uncoded SC-FDMA transmission
% SISO: simulation of several tx-rx user pairs without channel encoding
% Author: Uyen Ly DANG
% 30.07.2012: Basic SC-FDMA transmission chain
clear variables; close all;
clc;

%% System configuration
EbNo = [9 10]; % in dB

% User specific parameters
%           = [ user1 user2 user3 etc.]
User.bits_per_symb   = 2; %[2 2 2]; % Set number of bits per symbol: 2=QPSK, 4=16QAM=
User.R_c             = 1; %[1 1 1]; % Channel coder rate
User.M_mu            = 1; %[1 1 1];
User.N_nu            = 300; %[300 300 300];
User.var             = 1; %[1 1 1];
N_user = length(User.bits_per_symb); % Number of transmitter-receiver pairs

% SC-FDMA specific parameters
Scfdma.l_cp         = 144; % length of cyclic prefix
Scfdma.M            = 300; % length of DFT in transmitter
Scfdma.N            = 512; % length of IDFT in transmitter
Scfdma.nu_0         = 60;  % index of first occupied subcarrier (Localized mode)
Scfdma.N_scSymb     = 14;  % 2 slots with 7 SCFDMA-Symbols

% Channel parameters
Channel.type        =  'ITU-PB';'ITU-PA';
Channel.N_snapshot  = 1000;
Channel.N_sim       = 3; % Number of snapshots to be simulated
Channel.N_i         = 1; % Number of input streams per user(transmit antenna)
Channel.N_o         = 1; % Number of output streams per user(receive antenna)

Channel.start       = 1; 
Channel.offset      = [0 200 434] + Channel.start - 1; % starting point for

%% Sonstiges

% Get the postfix of the current file
current_postfix = getfilenamenumber(mfilename);

% Saving options
sim_file_name               = [ 'SC-FDMA_IA_' current_postfix ];
sim_file_name_not_finished  = [ 'Res_NF_' sim_file_name ] ;
sim_file_name_finished      = [ 'Res_F__' sim_file_name ] ;

% Initialize the random number generators with default values
randn('seed', 931316785) ; %#ok<RAND>     % for Matlab Version 
rand('seed', 931316785)  ;     %#ok<RAND>

% Get transmission channel
load([ 'LTE_channel_' Channel.type '_Anz' num2str(Channel.N_snapshot) '_cell_NR' num2str(Channel.N_o) '_NT' num2str(Channel.N_i) '.mat' ]);


% Get frame parameters, noise variance
idx_ch = zeros(N_user, Channel.N_snapshot);
num_bits = zeros(1,N_user); % N bits to be generated
EbNo_lin = 10.^(EbNo./10);
VarN = zeros(1,N_user);

for i_user = 1:N_user
    % channel index vector
    i_temp = 1:Channel.N_snapshot;
    idx_ch(i_user,:) = circshift(i_temp,[0 Channel.offset(i_user)]);
    
    % frame parameters
    num_bits(i_user)     =  Scfdma.M* User.bits_per_symb(i_user)*Scfdma.N_scSymb; %
end    

N = Scfdma.N;
M = Scfdma.M;
nu_0 = Scfdma.nu_0;
N_x = N+Scfdma.l_cp;

numSubframes = zeros(N_user, length(EbNo));

tic
for i_ebNo = 1:length(EbNo)
    % noise variance
    for i_user = 1:N_user
        VarN(i_user) = 1./(User.bits_per_symb(i_user).*User.R_c(i_user).*EbNo_lin(i_ebNo));
    end
    numSubframes(:, i_ebNo) = numSubframes(:, i_ebNo)+ones(i_user,1);

    for i_ch = 1:Channel.N_sim
        %% Transmitter
        % Signal generation
        tx_bits = cell(1,i_user);
        tx = cell(N_user,1);
        SymbTab = cell(1,N_user);
        for i_user = 1:N_user
            tx_bits{i_user} = randn(1,num_bits(i_user))>0;
            % Channel encoding
            % Channel interleaving
            % Bit mapping
            [tx{i_user}, SymbTab{i_user}] = bitMap(tx_bits{i_user}, User.bits_per_symb(i_user));
        end

        % SC-FDMA signal generation
        tx_sc    = cell(N_user,1);
        Tx       = zeros(1,M*Scfdma.N_scSymb);
        Tx_x     = zeros(1,N*Scfdma.N_scSymb);
        tx_x     = zeros(1,N*Scfdma.N_scSymb);
        tx_cp    = zeros(1,N*Scfdma.N_scSymb);
        
        % precoding vector, size (Nt, 1) where Nt is number of transmit
        % antennas        

        for i_user = 1:N_user
            tx_temp = tx{i_user};
            [V, G] = myPrecoding(h_cell, N_user, idx_ch, i_ch, VarN);
            for i_bl = 0:Scfdma.N_scSymb-1
                % Frequency domain: sc assignment and precoding
                Tx(i_bl*M+(1:M)) = fft(tx_temp(i_bl*M+(1:M)),M)./sqrt(M);        % V{i_user} * 
                Tx_x(i_bl*N+nu_0+(1:M)) = Tx(i_bl*M+(1:M));
                                
                % Time domain: insert cyclic prefix
                tx_x(i_bl*N+(1:N)) = ifft(Tx_x(i_bl*N+(1:N)),N).*(sqrt(N));
                tx_cp(i_bl*(N_x)+(1:(N_x))) = [tx_x(i_bl*N+(N-Scfdma.l_cp+1:N)) tx_x(i_bl*N+(1:N))];
            end
            tx_sc{i_user} = tx_cp;
        end
        clear Tx_temp; clear Txx_temp; clear txx_temp1; clear txx_temp2;

        %% Transmission
        rx_sc = cell(1,N_user);
        noise = cell(1,N_user);
        noise_var = zeros(1,N_user);
        
        for i_user = 1:N_user    
            rx_sc{i_user} = convLin(tx_sc{i_user},h_cell{idx_ch(i_user,i_ch)});
            noise{i_user} = sqrt(0.5*VarN(i_user)).*(randn(1,length(rx_sc{i_user}))+ 1i*randn(1,length(rx_sc{i_user})));
            %noise{i_user} = repmat(sqrt(0.5*noiseVar)',N_r,length(rx{i_user})).*(randn(N_r,length(rx{i_user}))+ 1i*randn(N_r,length(rx{i_user})));
            rx_sc{i_user} = rx_sc{i_user}+noise{i_user};
        end
  
        %% Receiver
        
        rx = cell(N_user,1);
        rx_bits = cell(N_user,1);
        Rx = zeros(1,M*Scfdma.N_scSymb);
        Rx_x = zeros(1,N*Scfdma.N_scSymb);
        Rx_e = zeros(1,N*Scfdma.N_scSymb);
        rx_x = zeros(1,N*Scfdma.N_scSymb);
        rx_cp = zeros(1,N*Scfdma.N_scSymb);
        for i_user = 1:N_user
            rx_temp = zeros(1,Scfdma.N_scSymb*M);
            %SC-FDMA Signal processing
            rx_cp = rx_sc{i_user};
            for i_bl = 0:Scfdma.N_scSymb-1                
                rx_x(i_bl*N+(1:N)) = rx_cp(i_bl*N_x+(Scfdma.l_cp+1:N_x));
                % Frequency domain: sc assignment/SC mapping
                Rx_x(i_bl*N+(1:N)) = fft(rx_x(i_bl*N+(1:N)),N)./sqrt(N);
                
                Rx(i_bl*M+(1:M)) = Rx_x(i_bl*N + nu_0+(1:M));
                
                % Equalization
                H  = fft(squeeze(h_cell{idx_ch(i_user,i_ch)}), M)./sqrt(M);
                Hh = H';
                F  = (Hh*H+ VarN(i_user))^-1 * Hh;
                
                Rx_e(i_bl*M+(1:M)) = F .* Rx(i_bl*M+(1:M));
                
                rx_temp(1,i_bl*M+(1:M)) = ifft(Rx_e(i_bl*M+(1:M)),M).*(sqrt(M));
                
            end
            rx{i_user} = rx_temp;
            % Channel Deinterleaving
            % Channel Encoding & Bit Demapping            
            [rx_bits{i_user},temp] = bitDemap(rx{i_user}, User.bits_per_symb(i_user), VarN(i_user));
        end
        % clear G; clear V;
        %% BLER/BER Calculation
        new_bitError = zeros(N_user,1);
        num_bit      = zeros(N_user,1);
        numBits      = zeros(N_user,length(EbNo));
        bitError     = zeros(N_user,length(EbNo));
        EbNo_sim     = zeros(1,length(EbNo));
        subframeError     = zeros(N_user,length(EbNo));
        % compute error rate
        for i_user = 1:N_user
            new_bitError(i_user)     = sum(tx_bits{i_user} ~= rx_bits{i_user} , 2);
            num_bit(i_user)          = length(rx_bits{i_user});
            numBits(i_user, i_ebNo)  = numBits(i_user, i_ebNo) + length(rx_bits{i_user});
        end

        new_subframeError        = new_bitError>0;
        bitError(:, i_ebNo)      = bitError(:, i_ebNo)+new_bitError;
        subframeError(:, i_ebNo) = subframeError(:, i_ebNo) + new_subframeError;
        EbNo_sim(i_ebNo)     = EbNo(i_ebNo);
        
        if ( i_ch == 1)
            est_finish_time = toc * Channel.N_snapshot ;
            days            = floor(est_finish_time/(3600*24));
            hours           = floor((est_finish_time-3600*24*days)/3600);
            minutes         = floor((est_finish_time-3600*24*days-3600*hours)/60);
            seconds         = round(est_finish_time-3600*24*days-3600*hours-60*minutes);
            fprintf('Estimated time untill finish of the simulation: %i d, %i h, %i m, %i s.\n', days, hours, minutes, seconds);
        end
        % Save not finished result after every 100 channel realizations
        if (mod(i_ch,100) == 0)
            BER     = bitError./numBits; %#ok<NASGU>
            BLER    = subframeError./numSubframes; %#ok<NASGU>
            rand_state  = rand('seed'); %#ok<RAND>
            randn_state = randn('seed'); %#ok<RAND>
            %save(sim_file_name_not_finished, 'i_ch', 'bitError', 'numBits', 'subframeError', 'numSubframes', 'BER','BLER','db_simulated','transmission_parameter','coder_input', 'sweap', 'N_user', 'rand_state', 'randn_state');
            save(sim_file_name_not_finished, 'i_ch', 'bitError', 'numBits', 'subframeError', 'numSubframes', 'BER','BLER','EbNo_sim', 'User','Channel','Scfdma', 'rand_state', 'randn_state');
        end
    end
    BER     = bitError./numBits;
    BLER    = subframeError./numSubframes;   
end

rand_state  = rand('seed'); %#ok<RAND>
randn_state = randn('seed'); %#ok<RAND>
%save(sim_file_name_finished,'subframeError','numSubframes','bitError','numBits','BER','BLER','EbNo_sim','transmission_parameter','coder_input', 'number_of_users','rand_state', 'randn_state');
save(sim_file_name_finished,'i_ch','subframeError','numSubframes','bitError','numBits','BER','BLER','EbNo_sim', 'User','Channel','Scfdma','rand_state', 'randn_state');
