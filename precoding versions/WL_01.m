%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Widely-Linear Transceiver Design for Interference Alignment                         
% Widely-linear filtering process the real and imaginary components of
% input symbols separately and independently.                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                   Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the list of input parameters which should be checked every time
% befor running the code.
% N_Users:                  total number of users
% User.bits_per_symb:       BPSK or QPSK...
% N_symbols:                number of symbols for one channel snapshot
% B:                        Number of parallel bits
% br:                       number of real-valued bits in B-bits
% nr:                       number of real-valued bits at the input of antenna
% Rs:                       augmented correlation matrix

% Note: For improper QPSK bitmapper and demapper should be modified
% properly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear variables; clear all; clc; % clear all variables and command window 

TEST = true;
%addpath('../MIMO_transmissionUpdated');
EbNo = [0:3:15]; % in dB

N_user = 3; % number of transmitter-receiver pairs

% User specific parameters
%                    = [user1 user2 user3 etc.]
User.bits_per_symb   = [2 2 2];%[2 2 2 2 2 2]; % Set number of bits per symbol: 2=QPSK, 4=16QAM=
User.R_c             = [1 1 1]; % Channel coder rate
User.var             = ones(1,N_user);

% Save File
current_postfix = getfilenamenumber(mfilename);
sim_file_name               = [ 'WL_BPSK_singlestream_' current_postfix];
sim_file_name_finished      = [ 'Res_F__' sim_file_name '_EbNo_' num2str(EbNo)] ;

% Channel parametersMIMO_03
Channel.N_snapshot  = 1000; %  Number of different snapshots per user
Channel.N_sim       = 1000; % Number of snapshots to be simulated
Channel.N_r         = 2; % Number of receive antennas
Channel.N_t         = 2; % Number of transmit antennas
Channel.N_i         = 1; % Number of input streams per user(T_tx)
Channel.N_o         = 3; % Number of output streams per user( N_rx * N_user)
Channel.test_channel = false;

if TEST
    Channel.N_sim = 100;
end

% Initialize the random number generators with default values
rng('shuffle');

% Get frame parameters, noise variance
idx_ch          = zeros(N_user, Channel.N_snapshot);
num_bits        = zeros(1,N_user); % N bits to be generated
N_symb          = 14; % must be multiple of N_tx;
EbNo_lin        = 10.^(EbNo./10);
varN            = zeros(1,N_user);

N_rx            = Channel.N_r;
N_tx            = Channel.N_t;

% Input Parameters
B               = 1; % number of symbols per streams 
br              = 0; % number of real valued symbols
Br              = 2*B - br;
nr              = 0; % real valued components at input of each antenna
Dr              = 2*N_tx - nr;

for i_user = 1:N_user   
    % number of bits per stream
    num_bits(i_user) = N_symb*User.bits_per_symb(i_user)*B;    
end

% Some initializations for BER/BLER evaluations
numSubframes    = zeros(N_user, length(EbNo));
bitError        = zeros(N_user,length(EbNo));
Sum_rate        = zeros(N_user,length(EbNo));
total_SumRate   = zeros(length(EbNo));
SINR            = zeros(N_user,length(EbNo));
new_bitError1   = zeros(N_user,length(EbNo));
EbNo_sim        = zeros(1,length(EbNo));
subframeError   = zeros(N_user,length(EbNo));
numBits         = zeros(N_user,length(EbNo));
total_power1    = zeros(N_user);
total_interferencePower     = zeros(N_user,1);


tic
for i_ebNo = 1:length(EbNo)
    new_bitError = zeros(N_user,1);
    
    for i_user = 1:N_user
        % noise variance 
        varN(i_user) = User.var(i_user)/(User.bits_per_symb(i_user).*User.R_c(i_user).*EbNo_lin(i_ebNo));
        
        % Variance of augmented symbol vector 
        if User.bits_per_symb(i_user) == 1 % When input symbols are real
            Rs(:,:,i_user) = eye(Br);
        elseif User.bits_per_symb(i_user) == 2 % When input symbols are complex
            Rs(:,:,i_user) = 0.5*eye(Br);
        end
    end
    numSubframes(:, i_ebNo) = numSubframes(:, i_ebNo)+ones(i_user,1);
    i_snapshot = 0;
    for i_ch = 1:Channel.N_sim 
       power = zeros(N_user,1); 
       
       % Transmission Channel
        h_all = zeros(N_rx,N_tx,N_user,N_user); % complex-valued
        H_all = zeros(2*N_rx,2*N_tx,N_user,N_user); % real-valued
        for i_tx = 1:N_user
            for i_rx = 1:N_user
                h_all(1:N_rx,1:N_tx,i_rx,i_tx) = (1/sqrt(N_tx))*[(randn(N_rx,N_tx) + 1i*randn(N_rx,N_tx))];
                % Real-valued transformation
                H_all(1:2*N_rx,1:2*N_tx,i_rx,i_tx) = [real(h_all(:,:,i_rx,i_tx)) -imag(h_all(:,:,i_rx,i_tx)) ; imag(h_all(:,:,i_rx,i_tx)) real(h_all(:,:,i_rx,i_tx))];
            end
        end

        
        %% Compute receive filter
        count               = 10;
        bfv                 = zeros(Dr,Br,N_user);
        bfv_New             = zeros(Dr,Br,N_user);
        bfv_norm            = zeros(count,1,N_user);
        rFilter_all         = zeros(2*N_rx,Br,N_user);
        rFilter_New         = zeros(2*N_rx,Br,N_user);
        innerSum_lemda      = zeros(N_rx,N_tx,N_user);
        lemda_New           = zeros(1,1,N_user);
        innerSum_MMSE_ini   = zeros(2*N_rx,2*N_rx,N_user);
        innerSum_MMSE_new   = zeros(N_rx,N_tx,N_user);
        lemda_check         = zeros(count,1,N_user);
        eta_sum            = zeros(100,1);
        
        % Initialization of beamforming vector 
        
        for i_tx = 1:N_user
            %bfv(:,:,i_tx) = randn(Dr,Br);
             [vector,value] = eig(h_all(:,:,i_tx,i_tx)'*h_all(:,:,i_tx,i_tx));
             if value(1,1) > value(2,2)
                 bfv1(:,:,i_tx) = vector(:,1);
            elseif value(1,1) < value(2,2)
                 bfv1(:,:,i_tx) = vector(:,2);
             end
             if Br == 1
               % bfv3(:,:,i_tx) = randn(Dr/2,Br) + 1i*randn(Dr/2,Br);
                bfv_ini(:,:,i_tx) =[real(bfv1(:,:,i_tx)) ; imag(bfv1(:,:,i_tx))];
             elseif Br > 1
               % bfv3(:,:,i_tx) = randn(Dr/2,Br/2) + 1i*randn(Dr/2,Br/2);
                bfv_ini(:,:,i_tx) =[real(bfv1(:,:,i_tx)) -imag(bfv1(:,:,i_tx)) ; imag(bfv1(:,:,i_tx)) real(bfv1(:,:,i_tx))];
            end
        end
        
        % Receive filter based on initial beamforming vector
        for i_rxr = 1:N_user
            for i_txt = 1:N_user
                 innerSum_MMSE_ini(:,:,i_rxr) = innerSum_MMSE_ini(:,:,i_rxr) + (H_all(:,:,i_rxr,i_txt)*bfv_ini(:,:,i_txt)*Rs(:,:,i_txt)*bfv_ini(:,:,i_txt)'*H_all(:,:,i_rxr,i_txt)');
            end
            rFilter_all(:,:,i_rxr) = (innerSum_MMSE_ini(:,:,i_rxr) + (varN(i_rxr)/2)*eye(2*N_rx))^(-1)*H_all(:,:,i_rxr,i_rxr)*bfv_ini(:,:,i_rxr)*Rs(:,:,i_rxr);
        end
        
        %% This loop is to update v_k and g_k
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here, beamforming vectors/matrices and receive filters are updated for
        % each user.
        % Note: beamformers and recieve filters are real-valued now
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:100 
            lemda_New = zeros(1,1,N_user);
            sum_temp = zeros(2*N_tx,2*N_tx,N_user);
            bfv_temp = zeros(Dr,Br,N_user);
            bfv_tempNorm = zeros(1,1,N_user);
            innerSum_lemda  = zeros(2*N_tx,2*N_tx,N_user);
            bfv_old(:,:,:) = bfv_New(:,:,:);
            
            %% First bfv_temp is computed, if trace(bfv_temp*Rs*bfv_temp^H) < = 1 then bfv_New = bfv_temp
            for i_tx1 = 1:N_user
                for i_rx1 = 1:N_user
                    sum_temp(:,:,i_tx1) = sum_temp(:,:,i_tx1) + H_all(:,:,i_rx1,i_tx1)'*rFilter_all(:,:,i_rx1)*Rs(:,:,i_rx1)*rFilter_all(:,:,i_rx1)'*H_all(:,:,i_rx1,i_tx1);
                end
                bfv_temp(:,:,i_tx1) = pinv(sum_temp(:,:,i_tx1))*H_all(:,:,i_tx1,i_tx1)'*rFilter_all(:,:,i_tx1)*Rs(:,:,i_tx1);
                bfv_tempNorm(:,:,i_tx1) = trace(bfv_temp(:,:,i_tx1)*Rs(:,:,i_tx1)*bfv_temp(:,:,i_tx1)');
            end
            
            % If trace(bfv_temp*Rs*bfv_temp^H) > 1 then bfv_New is computed here.
            for i_tx = 1:N_user
                if bfv_tempNorm(:,:,i_tx) <= 1
                    bfv_New(:,:,i_tx) = bfv_temp(:,:,i_tx);
                elseif bfv_tempNorm(:,:,i_tx) > 1
                    for i_rx = 1:N_user
                        innerSum_lemda(:,:,i_tx) = innerSum_lemda(:,:,i_tx) + H_all(:,:,i_rx,i_tx)'*rFilter_all(:,:,i_rx)*Rs(:,:,i_rx)*rFilter_all(:,:,i_rx)'*H_all(:,:,i_rx,i_tx);
                    end
                    
                    %Lambda Computation
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % This loop is for finding the lamda. Starting with 
                    % lemda_old =0, after each iteration lemda_old is replaced
                    % by lemda_New which is found in last iteration.
                    % beamforming vector is also computed for each lemda
                    % value. Loop will be terminated If Power constraint 
                    % (bfv_norm = 1) is fullfiled before maximum iteration 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    num1 = 0;
                    num = 0;
                    den = 0;
                    lemda_old = 0;
                    for i = 1:count % This loop is for finding lemda_k
                        num1 = pinv(innerSum_lemda(:,:,i_tx) + lemda_old*eye(2*N_tx)) * H_all(:,:,i_tx,i_tx)'*rFilter_all(:,:,i_tx)*Rs(:,:,i_tx);
                        num = trace(num1'*num1*Rs(:,:,i_tx))-1;
                        den = trace(2*Rs(:,:,i_tx)*num1'*pinv(innerSum_lemda(:,:,i_tx) + lemda_old*eye(2*N_tx))*num1);
                        lemda_New(:,:,i_tx) = lemda_old + ((num/den));
                        lemda_old = lemda_New(:,:,i_tx); % lemda_old is replace with lemda_New for next iteration
                        lemda_check(i,j,i_tx) = lemda_New(:,:,i_tx); % Should be positive (just to make sure that lemda values are right)
                        bfv_New(:,:,i_tx) = pinv(innerSum_lemda(:,:,i_tx) + (lemda_New(:,:,i_tx)*eye(2*N_tx))) * H_all(:,:,i_tx,i_tx)'*rFilter_all(:,:,i_tx)*Rs(:,:,i_tx);
                        bfv_norm(i,:,i_tx) = trace(bfv_New(:,:,i_tx)*Rs(:,:,i_tx)*bfv_New(:,:,i_tx)');
                        if bfv_norm(i,:,i_tx) == 1 % Power Constraint
                            break;
                        end
                    
                    end
                end
                norm_all(j,i_tx) = norm(bfv_New(:,:,i_tx),2)^2;
            end
            
            %% g_k_New
            % Recieve Filter Update (real-valued)
            innerSum_MMSE_new = zeros(2*N_rx,2*N_rx,N_user);
            Sum_MMSE1 = zeros(N_rx,N_tx,N_user);
            for i_rrx = 1:N_user
                for i_ttx = 1:N_user
                    innerSum_MMSE_new(:,:,i_rrx) = innerSum_MMSE_new(:,:,i_rrx) + H_all(:,:,i_rrx,i_ttx)*bfv_New(:,:,i_ttx)*Rs(:,:,i_ttx)*bfv_New(:,:,i_ttx)'*H_all(:,:,i_rrx,i_ttx)';
                end
                rFilter_New(:,:,i_rrx) = (innerSum_MMSE_new(:,:,i_rrx) + (varN(i_rrx)/2)*eye(2*N_rx))^(-1) * H_all(:,:,i_rrx,i_rrx)*bfv_New(:,:,i_rrx)*Rs(:,:,i_rrx);
                rFilter_all(:,:,i_rrx) = rFilter_New(:,:,i_rrx);
                %rFilter_normNew(:,:,i_rrx) = trace(rFilter_all(:,:,i_rrx)'*Rs(:,:,i_rrx)*rFilter_all(:,:,i_rrx));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                       Cross-Check
            % If this algorithem is working properly then sum MSE should be
            % decreasing with each iteration. eta should be decreasing if
            % it is not decreasing then algorithem is not working
            % properly. sum MSE is given here:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sum_MSE = zeros(Br,Br,N_user);
            eta = zeros(Br,Br);
            
            for ii_tx = 1:N_user
                for ii_rx = 1:N_user
                    if ii_tx ~= ii_rx
                        sum_MSE(:,:,ii_tx) = sum_MSE(:,:,ii_tx) + rFilter_New(:,:,ii_tx)'*H_all(:,:,ii_tx,ii_rx)*bfv_New(:,:,ii_rx);
                    end
                end
                eta(:,:) =(rFilter_New(:,:,ii_tx)'*H_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx)-eye(Br) + (sum_MSE(:,:,ii_tx)) + rFilter_New(:,:,ii_tx)'*rFilter_New(:,:,ii_tx)*(varN(ii_tx)/2)) *Rs(:,:,ii_tx)*... 
                                    (rFilter_New(:,:,ii_tx)'*H_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx)-eye(Br) + (sum_MSE(:,:,ii_tx)) + (rFilter_New(:,:,ii_tx)'*rFilter_New(:,:,ii_tx)*(varN(ii_tx)/2)))';
                eta_sum(j) = eta_sum(j)+ trace(eta(:,:));
            end
            
            % Convergence Limit
            Convergence_check(N_user,300) = 0;
            conver(300) = 0;
            for k = 1:N_user
                Convergence_check(k,j) = Convergence_check(k,j) + trace((bfv_New(:,:,k)-bfv_old(:,:,k))*(bfv_New(:,:,k)-bfv_old(:,:,k))');
                conver(j) = conver(j) + Convergence_check(k,j);
            end
            if conver(j) <= 10^-4
                break;
            end
        end
        
        %% Percentage of the Interference 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here, Percentage of Interference in the desired signal space is
        % computed. Formula used is: 
        % ro = (sum(lambda_j[Q_k]))/trace(Q_k), 
        % where lambda_j[Q_k] is the smallest eigen value of Q_k, and Q_k is the
        % interference (Q_k = sum_(j ~= k)(h(:,:,k,j)'bfv_new(:,j)bfv_new(:,j)'h(:,:,k,j))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        int_sum = zeros(N_rx,N_tx,N_user);
        int_sum1 = zeros(2*N_rx,2*N_rx,N_user);
        i=1;
            for i_rx = 1:N_user
                if i_rx == i+1
                    for i_tx = 1:i_rx
                        if i_rx ~= i_tx
                            int_sum1(:,:,i_rx) = int_sum1(:,:,i_rx) + H_all(:,:,i_rx,i_tx)*bfv_New(:,:,i_tx)*Rs(:,:,i_rx)*bfv_New(:,:,i_tx)'*H_all(:,:,i_rx,i_tx)';
                        end
                        
                    end
                    power1(i_rx) = abs((min(eig(int_sum1(:,:,i_rx)))+min(eig(int_sum1(:,:,i_rx))))/trace(int_sum1(:,:,i_rx)))*100;
                    total_power1(i_rx) = total_power1(i_rx) + power1(i_rx);
                    i = i+1;
                end
                
            end
            
            %%Sum Rate Computation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Here, Sum Rate achieved by all users is computed. Shannon
            % formula is used to compute rate:
            % R_k = 0.5*log_2(1 + SNR)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sum_Interference = zeros(Br,Br,N_user);
            for ii_tx = 1:N_user
                for ii_rx = 1:N_user
                    if ii_tx ~= ii_rx
                        sum_Interference(:,:,ii_tx) = sum_Interference(:,:,ii_tx) + Rs(:,:,ii_rx)*(abs(rFilter_New(:,:,ii_tx)'*H_all(:,:,ii_tx,ii_rx)*bfv_New(:,:,ii_rx))^2); % Interference  
                    end
                end
                desired_signal(:,:,ii_tx) = Rs(:,:,ii_tx)*(abs(rFilter_New(:,:,ii_tx)'*H_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx))^2); % desired signal
                Noisep(:,:,ii_tx) = (trace(rFilter_New(:,:,ii_tx)'*rFilter_New(:,:,ii_tx)))*((varN(ii_tx)/2)); % Noise
                R_k = (Br/2)*log2(1 + (trace(desired_signal(:,:,ii_tx))/(trace(sum_Interference(:,:,ii_tx))+Noisep(:,:,ii_tx)))); % Rate of individual user
                R_k(isnan(R_k)) = 0;
                total_SumRate(i_ebNo) = total_SumRate(i_ebNo) + R_k;
            end
        %% Transmitter
        % Signal generation
        % All users have the same signal length
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % On the transmitter side, first real and imaginary parts of the
        % input signal are separated. After real-valued transformation
        % transmit filter is applied which is also real-valued. After
        % applying transmit filter real and complex components are again
        % combine to transmit over the channel.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %if br == 0
            tx_bits             = zeros(1,num_bits(1),N_user);
            %tx                  = zeros(2*N_tx,N_symb,N_user);
            SymbTab             = cell(1,N_user);
            req_symbols         = zeros(1,N_symb,N_user);
            noise1               = zeros(N_rx,N_symb,N_user);
            noise               = zeros(2*N_rx,N_symb/B,N_user);
            rx                  = zeros(N_rx,N_symb/B,N_user);
            rx1                  = zeros(N_rx,N_symb,N_user);
            rx_req              = zeros(N_rx,N_symb,N_user);
            rx_int              = zeros(N_rx,N_symb,N_user);
        
            for i_user = 1:N_user            
                tx_bits(1,:,i_user) = randn(1,num_bits(i_user))>0;
            
                % Channel encoding
                % Channel interleaving
                % Bit mapping
                % bitMap takes only row vectors...          
                [tx_symbols(1,:,i_user),SymbTab{i_user}] = bitMap(tx_bits(1,:,i_user), User.bits_per_symb(i_user));
                if User.bits_per_symb(i_user) == 1
                    tx_symbols(1,:,i_user) = tx_symbols(1,:,i_user).*sqrt(User.var(i_user));
                    tx2(:,:,i_user) = reshape(tx_symbols(1,:,i_user) , [Br,N_symb/Br]); % arrange symbols to apply transmit filter
                    tx(:,:,i_user) = bfv_New(:,:,i_user)*tx2(:,:,i_user);
                elseif User.bits_per_symb(i_user) == 2
                    tx_symbols(1,:,i_user) = tx_symbols(1,:,i_user).*sqrt(User.var(i_user)); % Now the variance of each symbol is also half.
                    
                    % Antenna mapping: Input symbols are multiplied with
                    % beamforming vector/matrix
                    tx2(:,:,i_user) = reshape(tx_symbols(1,:,i_user) , [B,N_symb/B]);
                    tx2_aug(:,:,i_user) = [real(tx2(:,:,i_user)) ; imag(tx2(:,:,i_user))]; % Separation of real and imaginary part.
                    tx(:,:,i_user) = bfv_New(:,:,i_user)*tx2_aug(:,:,i_user); % Multiplication of augmented symbols with beamforming vector/matrix
                end
                for l = 1:N_symb % Combine again real and imaginary part to transmit over channel
                    for p =1:2
                        tx_com(p,l,i_user) = tx(p,l,i_user) + 1i*tx(p+2,l,i_user);
                    end
                end
            end
           for i_ruser = 1:N_user            
                % transmission of i_user to the receivers
                for i_tuser = 1:N_user
                    %rx1(:,:,i_ruser)= rx1(:,:,i_ruser) + convLin(tx1(:,:,i_tuser),h_all(:,:,i_ruser,i_tuser)); 
                    rx(:,:,i_ruser)= rx(:,:,i_ruser) + h_all(:,:,i_ruser,i_tuser)*tx_com(:,:,i_tuser); 
                end           
                noise1(:,:,i_ruser) = sqrt(0.5*(varN(i_ruser))).*(randn(N_rx,N_symb)+ 1i*randn(N_rx,N_symb));
                rx(:,:,i_ruser) = rx(:,:,i_ruser)+noise1(:,:,i_ruser);
            end

        %% Transmission
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % On the receiver side, first real and imaginary parts of the
        % received signal are separated. After real-valued transformation
        % received filter is applied which is also real-valued. After
        % applying receive filter real and complex components are again
        % combine to get the estimate.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           rx_symbols = zeros(N_symb,N_user);
           rx_bits     = zeros(1,N_symb*User.bits_per_symb(1),N_user);   
           for i_user = 1:N_user 
               % Antenna demapping
               rx_received(:,:,i_user) = [real(rx(:,:,i_user)) ; imag(rx(:,:,i_user))]; % Real and imaginary parts of received symbols are separated
               rx_augsymbols(:,:,i_user) = rFilter_New(:,:,i_user)'*rx_received(:,:,i_user); % receive filter is applied
               if br == 0 
                    for i = 1:N_symb/B % combine again real and imaginary parts to estimate the original symbols
                        for k = 1:B
                            rx_symbols1 (k,i,i_user) = (rx_augsymbols(k,i,i_user) + 1i*rx_augsymbols((Br/2)+k,i,i_user));
                        end
                    end
                   rx_symbols(:,i_user) = reshape(rx_symbols1(:,:,i_user),[N_symb,1]); % Reshaped for right order as on transmitter side
               elseif br ~= 0
                   % Antenna demapping
                   % rx_augsymbols(:,:,i_user) = rFilter_New(:,:,i_user)'*rx(:,:,i_user);
                   rx_symbols(:,i_user) = reshape(rx_augsymbols(:,:,i_user),[N_symb,1]);
               end
                
               % Symbol demapping
               [rx_bits(1,1:N_symb*User.bits_per_symb(i_user),i_user)] = HardDemap(rx_symbols(:,i_user),User.bits_per_symb(i_user));         
           end
           for i_user = 1:N_user
               new_bitError(i_user)     = sum(tx_bits(:,:,i_user) ~= rx_bits(:,:,i_user) , 2);            
               numBits(i_user, i_ebNo)  = numBits(i_user, i_ebNo) + length(rx_bits(:,:,i_user));
           end
           
           new_subframeError          = new_bitError>0;
           bitError(:, i_ebNo)      = bitError(:, i_ebNo)+new_bitError;
           subframeError(:, i_ebNo) = subframeError(:, i_ebNo) + new_subframeError;
           EbNo_sim(i_ebNo)     = EbNo(i_ebNo);
        
        if (mod(i_ch,10000) == 0)
            BER     = bitError./numBits; %#ok<NASGU>
            BLER    = subframeError./numSubframes; %#ok<NASGU>
            rand_state  = rand('seed'); %#ok<RAND>
            randn_state = randn('seed'); %#ok<RAND>
        end
    end
    BER(:,i_ebNo)     = bitError(:,i_ebNo)./numBits(:,i_ebNo);
    BLER(:,i_ebNo)    = subframeError(:,i_ebNo)./numSubframes(:,i_ebNo);
    total_interferencePower = total_power1./Channel.N_sim;
    total_sumRate(:,i_ebNo) = total_SumRate(i_ebNo)./Channel.N_sim;
   % check_sumRate(i_ebNo) = total_SumRate(i_ebNo)./Channel.N_sim;
    if ( i_ebNo == 1)
        est_finish_time = toc * length(EbNo) ;
        days            = floor(est_finish_time/(3600*24));
        hours           = floor((est_finish_time-3600*24*days)/3600);
        minutes         = floor((est_finish_time-3600*24*days-3600*hours)/60);
        seconds         = round(est_finish_time-3600*24*days-3600*hours-60*minutes);
        fprintf('Estimated time untill finish of the simulation: %i d, %i h, %i m, %i s.\n', days, hours, minutes, seconds);
    end
            
end


bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));
bit_err_theo = 0.5*erfc(sqrt(2*EbNo_lin)/sqrt(2));
figure;
%semilogy(EbNo, BER, 'g', EbNo, bit_err_theo, 'r', 'MarkerSize', 10, 'LineWidth', 2);
semilogy(EbNo, BER)
xlabel('Eb/No (dB)');
ylabel('Bit error rate');
title('Bit Error Rate');
% legend('Simulation','Theory');
grid on;
figure;
axis([1 N_user 0 35])
plot(1:1:N_user,total_interferencePower(1:N_user))
title('Percentage of interference in desired signal space')
xlabel('Number of Users');
ylabel('Percentage of Interference');
figure;
plot(EbNo,total_sumRate(1:length(EbNo)))
title('Sum Rate')
xlabel('EbNO (dB)');
ylabel('Sum Rate (bpcu)');
rand_state  = rand('seed'); %#ok<RAND>
randn_state = randn('seed'); %#ok<RAND>
save(sim_file_name_finished,'EbNo','check_sumRate','i_ch','subframeError','numSubframes','bitError','numBits','BER','BLER','EbNo_sim', 'User','Channel','rand_state', 'randn_state');