function [ V_new, G_new ] = myPrecodingWL( H_all, H_ch, VarN, N, WL )
%MYPRECODING_2 Summary of this function goes here
%   Detailed explanation goes here

%% Variables Initialization
N_user  = size(H_ch, 1);
Nr      = size(H_ch{1},1);
Nt      = size(H_ch{1},2);

G_new   = cell(N_user, 1);
G_all	= cell(N_user, 1);
v_init  = cell(N_user, 1);
V_init  = cell(N_user, 1);
V_new   = cell(N_user, 1);
V_old   = cell(N_user, 1);
V_norm  = cell(N_user, 1);
V_temp  = cell(N_user, 1);
V_tempnorm  = cell(N_user, 1);

sum_G_init  = cell(N_user, 1);
sum_G_new   = cell(N_user, 1);
sum_V       = cell(N_user, 1);
sumV_temp   = cell(N_user, 1);

sum_lambda  = cell(N_user, 1);
lambda_new  = cell(N_user, 1);
lambda_old  = cell(N_user, 1);

sum_MSE    = cell(N_user, 1);
epslon     = cell(N_user, 1);
% eta_sum    = cell(N_user, 1);

count = 10;
iterations = 1;

Convergence_check(iterations, N) = 0;
eta_sum = zeros(iterations, N);

Dr    = WL.Dr;
Br    = WL.Br;
Rs    = cell(N_user, 1);

%% Initializing v_MMSE_k and g_MMSE_k

for k_user = 1:N_user
    V_new{k_user}       = zeros(Dr, Br, N);
    G_new{k_user}       = zeros(2*Nt, Br, N);
    G_all{k_user}       = zeros(2*Nt, Br, N);
    sum_G_init{k_user}  = zeros(2*Nr, 2*Nt, N);
    sum_V{k_user}       = zeros(2*Nr, 2*Nt, N);
    V_tempnorm{k_user}  = zeros(1, N);
    V_norm{k_user}      = zeros(1, N);
end

%% Initializing Rs in frequency domain

for k_user = 1:N_user
    Rs{k_user}       = fft(WL.Rs(:,:,k_user), N, 3);
    %     eta_sum{k_user}  = zeros(1, iterations); % sum mean square error
end

%% Initializing beamforming vector
for idx = 1:N
    for k_user = 1:N_user
        [vector, value] = eig(H_ch{k_user, k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx));
        if value(1,1) > value(2,2)
            v_init{k_user}(:,:,idx) = vector(:,1);
        elseif value(1,1) < value(2,2)
            v_init{k_user}(:,:,idx) = vector(:,2);
        end
        
        if Br == 1
            V_init{k_user}(:,:,idx) =[real(v_init{k_user}(:,:,idx)) ; imag(v_init{k_user}(:,:,idx))];
        elseif Br > 1
            V_init{k_user}(:,:,idx) =[real(v_init{k_user}(:,:,idx)), -imag(v_init{k_user}(:,:,idx)) ; imag(v_init{k_user}(:,:,idx)), real(v_init{k_user}(:,:,idx))];
        end
    end
end

%% update v_MMSE_k and g_MMSE_k
for idx = 1:N
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_G_init{k_user}(:,:,idx)  = sum_G_init{k_user}(:,:,idx) + H_all{k_user, j_user}(:,:,idx)* V_init{j_user}(:,:,idx)* Rs{j_user}(:,:,idx) * V_init{j_user}(:,:,idx)' * H_all{k_user, j_user}(:,:,idx)';
        end
        G_all{k_user}(:,:,idx) = pinv(sum_G_init{k_user}(:,:,idx) + (VarN(k_user)/2)*eye(2*Nr)) * H_all{k_user, k_user}(:,:,idx) * V_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
    end
end


for idx = 1:N
    for j = 1:iterations
        for k_user = 1:N_user
            lambda_new{k_user}      = zeros(1,N);
            lambda_old{k_user}      = zeros(1,N);
            sumV_temp{k_user}       = zeros(Dr, Dr, N);
            V_temp{k_user}          = zeros(Dr, Br, N);
            sum_lambda{k_user}      = zeros(Dr, Dr, N);
            
            sum_MSE{k_user} = zeros(Br, Br, N);
            %         epslon{k_user}  = zeros(j, N);
        end
        
        for k_user = 1:N_user
            V_old{k_user}(:,:,idx)  = V_new{k_user}(:,:,idx);
            for j_user = 1:N_user
                sumV_temp{k_user}(:,:,idx)  = sumV_temp{k_user}(:,:,idx) + H_all{j_user, k_user}(:,:,idx)' * G_all{j_user}(:,:,idx) * Rs{j_user}(:,:,idx) * G_all{j_user}(:,:,idx)' * H_all{j_user, k_user}(:,:,idx);
            end
            V_temp{k_user}(:,:,idx)   = pinv(sumV_temp{k_user}(:,:,idx)) * H_all{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
            V_tempnorm{k_user}(1,idx) = trace(V_temp{k_user}(:,:,idx) * Rs{k_user}(:,:,idx) * V_temp{k_user}(:,:,idx)');
        end
        
        for k_user = 1:N_user
            if V_tempnorm{k_user}(1,idx) <= 1
                V_new{k_user}(:,:,idx) = V_init{k_user}(:,:,idx);
            elseif V_tempnorm{k_user}(1,idx) > 1
                for j_user = 1:N_user
                    sum_lambda{k_user}(:,:,idx)  = sum_lambda{k_user}(:,:,idx) + H_all{j_user, k_user}(:,:,idx)' * G_all{j_user}(:,:,idx) * Rs{j_user}(:,:,idx) * G_all{j_user}(:,:,idx)' * H_all{j_user, k_user}(:,:,idx);
                end
                
                lambda_old{k_user}(1,idx) =  0;
                
                for i = 1:count
                    inverse_mat   = pinv(sum_lambda{k_user}(:,:,idx) + lambda_old{k_user}(1,idx)*eye(2*Nt));
                    inner_sum_num = inverse_mat * inverse_mat;
                    inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                    
                    num = trace( Rs{k_user}(:,:,idx) * G_all{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx)* inner_sum_num * H_all{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * Rs{k_user}(:,:,idx)) - 1;
                    den = trace( 2 * Rs{k_user}(:,:,idx) * G_all{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx) * inner_sum_den * H_all{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * Rs{k_user}(:,:,idx));
                    
                    lambda_new{k_user}(1,idx)   = lambda_old{k_user}(1,idx) + (num/den);
                    lambda_old{k_user}(1,idx)   = lambda_new{k_user}(1,idx);
                    
                    V_new{k_user}(:,:,idx)      = pinv(sum_lambda{k_user}(:,:,idx) + lambda_new{k_user}(1,idx)*eye(2*Nt)) * H_all{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
                    V_norm{k_user}(1,idx)       = trace(V_new{k_user}(:,:,idx) *  Rs{k_user}(:,:,idx) * V_new{k_user}(:,:,idx)');
                    
                    if V_norm{k_user}(1,idx) == 1
                        break;
                    end
                end
            end
        end
        
        for k_user = 1:N_user
            sum_G_new{k_user}(:,:,idx)  = zeros(2*Nr, 2*Nt);
            for j_user = 1:N_user
                sum_G_new{k_user}(:,:,idx)  = sum_G_new{k_user}(:,:,idx) + H_all{k_user, j_user}(:,:,idx) * V_new{j_user}(:,:,idx) * Rs{j_user}(:,:,idx) * V_new{j_user}(:,:,idx)'*H_all{k_user, j_user}(:,:,idx)';
            end
            G_new{k_user}(:,:,idx) = pinv(sum_G_new{k_user}(:,:,idx) + (VarN(k_user)/2)*eye(2*Nr)) * H_all{k_user, k_user}(:,:,idx) * V_new{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
            G_all{k_user}(:,:,idx) = G_new{k_user}(:,:,idx);
        end
        
        for k_user = 1:N_user
            for j_user = 1:N_user
                if k_user ~= j_user
                    sum_MSE{k_user}(:,:,idx) = sum_MSE{k_user}(:,:,idx) + G_all{k_user}(:,:,idx)' * H_all{k_user, j_user}(:,:,idx) * V_new{j_user}(:,:,idx);
                end
            end
            epslon{k_user}(:,:) = (((G_all{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx) * V_new{k_user}(:,:,idx)) - eye(Br)) + sum_MSE{k_user}(:,:,idx) + G_all{k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * (VarN(k_user)/2)) *...
                Rs{k_user}(:,:,idx) * ...
                (((G_all{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx) * V_new{k_user}(:,:,idx)) - eye(Br)) + sum_MSE{k_user}(:,:,idx) + G_all{k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) * (VarN(k_user)/2))';
            
            eta_sum(j, idx) = eta_sum(j, idx) + trace(epslon{k_user}(:,:));
        end
        
        
        for k_user = 1:N_user
            Convergence_check(j, idx) = Convergence_check(j, idx) + trace( (V_new{k_user}(:,:,idx) - V_old{k_user}(:,:,idx)) * (V_new{k_user}(:,:,idx) - V_old{k_user}(:,:,idx))' );
        end
        if Convergence_check(j, idx) <= 10^-4
            break;
        end
    end
end
end