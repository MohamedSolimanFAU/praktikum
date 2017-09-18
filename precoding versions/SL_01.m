%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Strictly Linear Transceiver Design for Interference Alignment        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                   Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here is the list of input parameters which should be checked every time
% befor running the code.
% N_Users:                  total number of users
% User.bits_per_symb:       BPSK or QPSK...
% N_symbols:                number of symbols for one channel snapshot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear variables; clear all; clc; % clear all variables and command window 

TEST = true;
%addpath('../MIMO_transmissionUpdated');

EbNo = [0:3:15]; % in dB

N_user = 3; % number of transmitter-receiver pairs

% User specific parameters
%                    = [user1 user2 user3 etc.]
User.bits_per_symb   = [2 2 2]; % Set number of bits per symbol: 2=QPSK, 4=16QAM=
User.R_c             = [1 1 1]; % Channel coder rate
User.var             = ones(1,N_user);

%save file
current_postfix             = getfilenamenumber(mfilename);
sim_file_name               = [ 'SL_QPSK_Singlestream_' current_postfix ];
sim_file_name_finished      = [ 'Res_F__' sim_file_name '_EbNo_' num2str(EbNo)] ;

% Channel parametersMIMO_03
Channel.N_snapshot  = 1000; %  Number of different snapshots per user
Channel.N_sim       = 10000; % Number of snapshots to be simulated
Channel.N_r         = 1; % Number of receive antennas
Channel.N_t         = 1; % Number of transmit antennas
Channel.N_i         = 1; % Number of input streams per user(T_tx)
Channel.N_o         = 3; % Number of output streams per user( N_rx * N_user)
Channel.test_channel = false;

if TEST
    Channel.N_sim = 100;
end

% Initialize the random number generators with default values
rng default;

% Get frame parameters, noise variance
idx_ch          = zeros(N_user, Channel.N_snapshot);
num_bits        = zeros(1,N_user); % N bits to be generated
N_symb          = 14; % must be multiple of N_tx;
EbNo_lin        = 10.^(EbNo./10);
varN            = zeros(1,N_user);

N_rx            = Channel.N_r;
N_tx            = Channel.N_t;

for i_user = 1:N_user   
    % number of bits per stream
    num_bits(i_user) = N_symb*User.bits_per_symb(i_user);    
end

% Some initializations for BER/BLER evaluations
numSubframes    = zeros(N_user, length(EbNo));
bitError        = zeros(N_user,length(EbNo));
Sum_rate        = zeros(N_user,length(EbNo));
new_bitError1   = zeros(N_user,length(EbNo));
EbNo_sim        = zeros(1,length(EbNo));
subframeError   = zeros(N_user,length(EbNo));
numBits         = zeros(N_user,length(EbNo));
total_interferencePower = zeros(N_user,1);
total_power1 = zeros(N_user,1);
tic
for i_ebNo = 1:length(EbNo)
    new_bitError = zeros(N_user,1);
    
    for i_user = 1:N_user
        % noise variance 
        varN(i_user) = User.var(i_user)/(User.bits_per_symb(i_user).*User.R_c(i_user).*EbNo_lin(i_ebNo));
    end
    numSubframes(:, i_ebNo) = numSubframes(:, i_ebNo)+ones(i_user,1);
    i_snapshot = 0;
    for i_ch = 1:Channel.N_sim 
       
       % Transmission Channel
        h_all = zeros(N_rx,N_tx,N_user,N_user);  
        for i_tx = 1:N_user
            for i_rx = 1:N_user
                h_all(1:N_rx,1:N_tx,i_rx,i_tx) = 1/sqrt(N_tx)*[randn(N_rx,N_tx) + 1i*randn(N_rx,N_tx)];
            end
        end

        
        %% Computation of Beamforming vecotors and Receive filters
        %initialization of variables requried for tranceivers
        count               = 10;
        bfv_ini             = zeros(N_rx,1,N_user); % Beamforming Vector
        bfv_New             = zeros(N_rx,1,N_user);
        rFilter_all         = zeros(N_rx,1,N_user); % Receive Filter
        rFilter_New         = zeros(N_rx,1,N_user);
        lemda_old           = zeros(1,1,N_user);
        lemda_New           = zeros(1,1,N_user);
        lemda_check         = zeros(count,1,N_user);
        innerSum_lemda      = zeros(N_rx,N_tx,N_user);
        innerSum_MMSE_ini   = zeros(N_rx,N_tx,N_user);
        innerSum_MMSE_new   = zeros(N_rx,N_tx,N_user);
        bfv_norm            = zeros(count,1,N_user);
        
        
        % Initialization of beamforming vector 
        for i_tx = 1:N_user
            [vector,value] = eig(h_all(:,:,i_tx,i_tx)'*h_all(:,:,i_tx,i_tx));
%             if value(1,1) > value(2,2)
%                 bfv_ini(:,:,i_tx) = vector(:,1);
%             elseif value(1,1) < value(2,2)
%                 bfv_ini(:,:,i_tx) = vector(:,2);
%             end
            bfv_ini(:,:,i_tx) = randn(N_rx,1) + 1i*randn(N_rx,1);
        end
       
        
        % Receive filter based on initial beamforming vector
        for i_rxr = 1:N_user
            for i_txt = 1:N_user
                 innerSum_MMSE_ini(:,:,i_rxr) = innerSum_MMSE_ini(:,:,i_rxr) + (h_all(:,:,i_rxr,i_txt)*bfv_ini(:,:,i_txt)*bfv_ini(:,:,i_txt)'*h_all(:,:,i_rxr,i_txt)');
            end
            rFilter_all(:,:,i_rxr) = (innerSum_MMSE_ini(:,:,i_rxr) + (varN(i_rxr)*eye(N_rx)))^(-1)*h_all(:,:,i_rxr,i_rxr)*bfv_ini(:,:,i_rxr);
            rFilter_norm(:,:,i_rxr) = norm(rFilter_all(:,:,i_rxr),2)^2;
        end
        % check whether beamforing vectors and receive filters are
        % initialize properly or not. checki should be close to one.
        for i_check = 1:N_user
            checki(:,:,i_check) = rFilter_all(:,:,i_check)'*h_all(:,:,i_check,i_check)*bfv_ini(:,:,i_check);
        end
        
        eta_sum = zeros(100,1); % sum mean square error
        
        %% This loop is to update v_k and g_k
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here, beamforming vectors and receive filters are updated for
        % each user
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:100
            lemda_old = zeros(1,1,N_user);
            lemda_New = zeros(1,1,N_user);
            sum_temp = zeros(N_tx,N_rx,N_user);
            bfv_temp = zeros(N_tx,1,N_user);
            bfv_tempNorm = zeros(1,1,N_user);
            innerSum_lemda  = zeros(N_rx,N_tx,N_user);
            bfv_old(:,:,:) = bfv_New(:,:,:);
            
            %% First bfv_temp is computed, if ||bfv_temp||^2 < = 1 then bfv_New = bfv_temp
            for i_tx1 = 1:N_user
                for i_rx1 = 1:N_user
                    sum_temp(:,:,i_tx1) = sum_temp(:,:,i_tx1) + h_all(:,:,i_rx1,i_tx1)'*rFilter_all(:,:,i_rx1)*rFilter_all(:,:,i_rx1)'*h_all(:,:,i_rx1,i_tx1);
                end
                bfv_temp(:,:,i_tx1) = pinv(sum_temp(:,:,i_tx1))*h_all(:,:,i_tx1,i_tx1)'*rFilter_all(:,:,i_tx1);
                bfv_tempNorm(:,:,i_tx1) = norm(bfv_temp(:,:,i_tx1),2)^2;
            end
            
            % If ||bfv_temp||^2 > 1 then bfv_New is computed here.
            for i_tx = 1:N_user
                if bfv_tempNorm(:,:,i_tx) <= 1
                    bfv_New(:,:,i_tx) = bfv_temp(:,:,i_tx);
                elseif bfv_tempNorm(:,:,i_tx) > 1
                    for i_rx = 1:N_user
                        innerSum_lemda(:,:,i_tx) = innerSum_lemda(:,:,i_tx) + h_all(:,:,i_rx,i_tx)'*rFilter_all(:,:,i_rx)*rFilter_all(:,:,i_rx)'*h_all(:,:,i_rx,i_tx);
                    end
                    num = 0;
                    den = 0;
                    lemda_old = 0;
                    
                    %Lambda Computation
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % This loop is for finding the lamda. Starting with 
                    % lemda_old =0, after each iteration lemda_old is replaced
                    % by lemda_New which is found in last iteration.
                    % beamforming vector is also computed for each lemda
                    % value. Loop will be terminated If Power constraint 
                    % (bfv_norm = 1) is fullfiled before maximum iteration 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    for i = 1:count 
                        num = rFilter_all(:,:,i_tx)'*h_all(:,:,i_tx,i_tx)*((innerSum_lemda(:,:,i_tx) + (lemda_old*eye(N_rx)))^(-2))*h_all(:,:,i_tx,i_tx)'*rFilter_all(:,:,i_tx)-1;
                        den = 2*rFilter_all(:,:,i_tx)'*h_all(:,:,i_tx,i_tx)*((innerSum_lemda(:,:,i_tx) + (lemda_old*eye(N_rx)))^(-3))*h_all(:,:,i_tx,i_tx)'*rFilter_all(:,:,i_tx);
                        lemda_New(:,:,i_tx) = lemda_old + (num/den);
                        lemda_old = lemda_New(:,:,i_tx); % lemda_old is replace with lemda_New
                        lemda_check(i,j,i_tx) = lemda_New(:,:,i_tx); % lemda_check should be positive
                        bfv_New(:,:,i_tx) = (innerSum_lemda(:,:,i_tx) + (lemda_New(:,:,i_tx)*eye(N_rx)))^(-1) * h_all(:,:,i_tx,i_tx)'*rFilter_all(:,:,i_tx);
                        bfv_norm(i,:,i_tx) = norm(bfv_New(:,:,i_tx),2)^2;
                        if bfv_norm(i,:,i_tx) == 1 % Power Constraint 
                            break;
                        end
                    
                    end
                end
                norm_all(j,i_tx) = norm(bfv_New(:,:,i_tx),2)^2;
            end
            
            %% g_k_New
            % Recieve Filter Update
            innerSum_MMSE_new = zeros(N_rx,N_tx,N_user);
            for i_rrx = 1:N_user
                for i_ttx = 1:N_user
                    innerSum_MMSE_new(:,:,i_rrx) = innerSum_MMSE_new(:,:,i_rrx) + h_all(:,:,i_rrx,i_ttx)*bfv_New(:,:,i_ttx)*bfv_New(:,:,i_ttx)'*h_all(:,:,i_rrx,i_ttx)';
                end
                rFilter_New(:,:,i_rrx) = (innerSum_MMSE_new(:,:,i_rrx) + (varN(i_rrx)*eye(N_rx)))^(-1) * h_all(:,:,i_rrx,i_rrx)*bfv_New(:,:,i_rrx);
                rFilter_all(:,:,i_rrx) = rFilter_New(:,:,i_rrx);
                rFilter_normNew(:,:,i_rrx) = norm(rFilter_all(:,:,i_rrx),2)^2;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If this algorithem is working properly then sum MSE should be
            % decreasing with each iteration. 
            % sum MSE is given here:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sum_MSE = zeros(1,1,N_user);
            eta = zeros(j,N_user); % mean square error
            
            for ii_tx = 1:N_user
                for ii_rx = 1:N_user
                    if ii_tx ~= ii_rx
                        sum_MSE(:,:,ii_tx) = sum_MSE(:,:,ii_tx) + abs(rFilter_New(:,:,ii_tx)'*h_all(:,:,ii_tx,ii_rx)*bfv_New(:,:,ii_rx))^2;
                    end
                end
                eta(j,ii_tx) = abs(rFilter_New(:,:,ii_tx)'*h_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx)-1)^2 + sum_MSE(:,:,ii_tx) + (norm(rFilter_New(:,:,ii_tx),2)^2)*varN(ii_tx);
                eta_sum(j) = eta_sum(j)+ eta(j,ii_tx); % Sum mean square error
            end
            % Convergence limit: 
            Convergence_check(300) = 0;
            for k = 1:N_user
                Convergence_check(j) = Convergence_check(j) + norm(bfv_New(:,:,k)-bfv_old(:,:,k),2);
            end
            if Convergence_check(j) <= 10^-4
                break;
            end
        end
        for i_check = 1:N_user
            check_update(:,:,i_check) = rFilter_New(:,:,i_check)'*h_all(:,:,i_check,i_check)*bfv_New(:,:,i_check);
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
        int_sum1 = zeros(N_rx,N_tx,N_user);
        i=1;
            for i_rx = 1:N_user
                if i_rx == i+1
                    for i_tx = 1:i_rx
                        if i_rx ~= i_tx
                            int_sum1(:,:,i_rx) = int_sum1(:,:,i_rx) + h_all(:,:,i_rx,i_tx)*bfv_New(:,:,i_tx)*bfv_New(:,:,i_tx)'*h_all(:,:,i_rx,i_tx)';
                        end
                        power1(i_rx) = min(eig(int_sum1(:,:,i_rx)))/trace(int_sum1(:,:,i_rx))*100;
                        total_power1(i_rx) = total_power1(i_rx) + power1(i_rx); % Interference power of each user
                    end
                    power1(i_rx) = min(eig(int_sum1(:,:,i_rx)))/trace(int_sum1(:,:,i_rx))*100;
                    total_power1(i_rx) = total_power1(i_rx) + power1(i_rx); 
                    i = i+1;
                end   
            end
            
            %%Sum Rate Computation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Here, Sum Rate achieved by all users is computed. Shannon
            % formula is used to compute rate:
            % R_k = log_2(1 + SNR) for complex
            % R_k = log_2(1 + SNR) for real
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sum_Interference = zeros(1,1,N_user);
            for ii_tx = 1:N_user
                for ii_rx = 1:N_user
                    if ii_tx ~= ii_rx
                        sum_Interference(:,:,ii_tx) = sum_Interference(:,:,ii_tx) + abs(rFilter_New(:,:,ii_tx)'*h_all(:,:,ii_tx,ii_rx)*bfv_New(:,:,ii_rx))^2;
                    end
                end
                desired_signal = abs(rFilter_New(:,:,ii_tx)'*h_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx))^2;
                interference_Noise = sum_Interference(:,:,ii_tx) + norm(rFilter_New(:,:,ii_tx),2)^2*varN(ii_tx);
                Sum_rate(ii_tx,i_ebNo) = Sum_rate(ii_tx,i_ebNo) + 0.5*log2(1+((abs(rFilter_New(:,:,ii_tx)'*h_all(:,:,ii_tx,ii_tx)*bfv_New(:,:,ii_tx))^2)/(sum_Interference(:,:,ii_tx) + norm(rFilter_New(:,:,ii_tx),2)^2*(varN(ii_tx)))));
            end

        %% Transmitter
        % Signal generation
        % All users have the same signal length
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % On the transmitter side, transmit filter is applied on input
        % symbos. After applying transmit filter the signal transmitted
        % over channel.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tx_bits             = zeros(1,num_bits(1),N_user);
        tx                  = zeros(N_tx,N_symb,N_user);
        SymbTab             = cell(1,N_user);
        req_symbols         = zeros(1,N_symb,N_user);
        int_symbols         = zeros(1,N_symb,N_user);
        noise               = zeros(N_rx,N_symb,N_user);
        rx                  = zeros(N_rx,N_symb,N_user);
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
            tx_symbols(1,:,i_user) = tx_symbols(1,:,i_user)*sqrt(User.var(i_user));
            
            % Antenna mapping
            for p =1:N_symb
                tx(:,p,i_user) = bfv_New(:,:,i_user)*tx_symbols(1,p,i_user);
            end
   
        end
       
        for i_ruser = 1:N_user            
            % transmission of i_user to the receivers
            for i_tuser = 1:N_user
                %rx(:,:,i_ruser)= rx(:,:,i_ruser) + convLin(tx(:,:,i_tuser),h_all(:,:,i_ruser,i_tuser)); 
                rx(:,:,i_ruser)= rx(:,:,i_ruser) + h_all(:,:,i_ruser,i_tuser)*tx(:,:,i_tuser); 
            end           
            noise(:,:,i_ruser) = sqrt(0.5*varN(i_ruser)).*(randn(N_rx,N_symb)+ 1i*randn(N_rx,N_symb)); 
            rx(:,:,i_ruser) = rx(:,:,i_ruser)+noise(:,:,i_ruser);
        end

        %% Transmission
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % On the receiver side, receive filter is applied to estimate the
        % symbols
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           rx_symbols = zeros(N_symb,N_user);
           rx_bits     = zeros(1,N_symb*User.bits_per_symb(1),N_user);
           
           for i_user = 1:N_user 
               % Antenna demapping: received signal is multiplied with
               % receive filter
               rx_symbols(:,i_user) = rFilter_New(:,:,i_user)'*rx(:,:,i_user);
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
    total_sumRate(:,i_ebNo) = Sum_rate(:,i_ebNo)./Channel.N_sim;
%     Sum_rate(:,i_ebNo) = 
    if ( i_ebNo == 1)
        est_finish_time = toc * length(EbNo) ;
        days            = floor(est_finish_time/(3600*24));
        hours           = floor((est_finish_time-3600*24*days)/3600);
        minutes         = floor((est_finish_time-3600*24*days-3600*hours)/60);
        seconds         = round(est_finish_time-3600*24*days-3600*hours-60*minutes);
        fprintf('Simulation Started at: %i.%i.%i ; %i:%i:%i.\n', fix(clock));
        fprintf('Estimated time untill finish of the simulation: %i d, %i h, %i m, %i s.\n', days, hours, minutes, seconds);
    end
            
end


bit_err_theo = 0.5*erfc(sqrt(EbNo_lin));
figure;
%semilogy(EbNo, BER, 'g', EbNo, bit_err_theo, 'r', 'MarkerSize', 10, 'LineWidth', 2);
semilogy(EbNo, BER)
xlabel('Eb/No (dB)');
ylabel('Bit error rate');
title('BPSK bit error rate');
% legend('Simulation','Theory');
grid on;
figure;
axis([1 N_user 0 35])
plot(1:1:N_user,total_interferencePower(1:N_user))
title(['N-users = ' ,N_user, 'N-tx = N-rx = ' ,N_tx])
xlabel('Number of Users');
ylabel('Percentage of Interference');
figure;
plot(EbNo,total_sumRate(:,1:length(EbNo)))
%title('')
xlabel('EbNO (dB)');
ylabel('Rate (bpcu)');
rand_state  = rand('seed'); %#ok<RAND>
randn_state = randn('seed'); %#ok<RAND>
save(sim_file_name_finished,'i_ch','subframeError','numSubframes','bitError','numBits','BER','BLER','EbNo_sim', 'User','Channel','rand_state', 'randn_state');