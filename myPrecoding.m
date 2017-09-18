function [ V_new, G ] = myPrecoding( H_ch, VarN, N )
%MYPRECODING_2 Summary of this function goes here
%   Detailed explanation goes here

%% Variables Initialization
N_user  = size(H_ch, 1);
Nr      = size(H_ch{1},1);
Nt      = size(H_ch{1},2);

G       = cell(N_user, 1);
G_all	= cell(N_user, 1);
V_new   = cell(N_user, 1);
V_old   = cell(N_user, 1);
v_init  = cell(N_user, 1);
V_init  = cell(N_user, 1);
V_temp  = cell(N_user, 1);
V_norm  = cell(N_user, 1);
V_tempnorm  = cell(N_user, 1);

sum_G      = cell(N_user, 1);
sum_V      = cell(N_user, 1);
sumV_temp  = cell(N_user, 1);
sum_lambda = cell(N_user, 1);
lambda_new = cell(N_user, 1);
lambda_old = cell(N_user, 1);

sum_MSE    = cell(N_user, 1);
epslon     = cell(N_user, 1);

count = 10;
iterations = 1000;
%% Initializing v_MMSE_k & g_MMSE_k

for k_user = 1:N_user
    v_init{k_user}     = randn(Nt, 1)+ 1i* randn(Nt, 1);
    V_init{k_user}     = fft(v_init{k_user}, N, 3);
    V_new{k_user}      = zeros(Nt, 1, N);
    G{k_user}          = zeros(Nt, 1, N);
    G_all{k_user}      = zeros(Nt, 1, N);
    sum_G{k_user}      = zeros(Nr, Nt, N);
    sum_V{k_user}      = zeros(Nr, Nt, N);
    V_tempnorm{k_user} = zeros(Nr, Nt, N);
end

Convergence_check(iterations, N) = 0;

for idx = 1:N
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_ch{k_user, j_user}(:,:,idx)*V_init{j_user}(:,:,idx)*V_init{j_user}(:,:,idx)'*H_ch{k_user, j_user}(:,:,idx)';
        end
        G_all{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + VarN(k_user)*eye(Nr)) * H_ch{k_user, k_user}(:,:,idx) * V_init{k_user}(:,:,idx);
    end
end

eta_sum = zeros(1, iterations); % sum mean square error

for j = 1:iterations
    for k_user = 1:N_user
        lambda_new{k_user}      = zeros(1,N);
        lambda_old{k_user}      = zeros(1,N);
        sumV_temp{k_user}       = zeros(Nr, Nt, N);
        V_temp{k_user}          = zeros(Nt, 1, N);
        sum_lambda{k_user}      = zeros(Nr, Nt, N);
    end
    
    for k_user = 1:N_user
        sum_MSE{k_user} = zeros(1, 1, N);
        epslon{k_user}  = zeros(j, N);
    end
    
    for idx = 1:N
        
        
        for k_user = 1:N_user
            
            V_old{k_user}(:,:,idx)  = V_new{k_user}(:,:,idx);
            
            for j_user = 1:N_user
                sumV_temp{k_user}(:,:,idx)  = sumV_temp{k_user}(:,:,idx) + H_ch{j_user, k_user}(:,:,idx)'*G_all{j_user}(:,:,idx)*G_all{j_user}(:,:,idx)'*H_ch{j_user, k_user}(:,:,idx);
            end
            V_temp{k_user}(:,:,idx) = pinv(sumV_temp{k_user}(:,:,idx)) * H_ch{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx);
            V_tempnorm{k_user}(:,:,idx) = norm(V_temp{k_user}(:,:,idx), 2)^2;
        end
        
        for k_user = 1:N_user
            num = 0;
            den = 0;
            lambda_old{k_user}(1,idx) =  0;
            if V_tempnorm{k_user}(:,:,idx) <= 1
                V_new{k_user}(:,:,idx) = V_temp{k_user}(:,:,idx);
            elseif V_tempnorm{k_user}(:,:,idx) > 1
                for j_user = 1:N_user
                    sum_lambda{k_user}(:,:,idx)  = sum_lambda{k_user}(:,:,idx) + H_ch{j_user, k_user}(:,:,idx)'*G_all{j_user}(:,:,idx)*G_all{j_user}(:,:,idx)'*H_ch{j_user, k_user}(:,:,idx);
                end
                
                
                for i = 1:count
                    inverse_mat   = pinv(sum_lambda{k_user}(:,:,idx) + lambda_new{k_user}(1,idx)*eye(Nt));
                    inner_sum_num = inverse_mat * inverse_mat;
                    inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                    
                    num =     G_all{k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx)* inner_sum_num * H_ch{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx) - 1;
                    den = 2 * G_all{k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx)* inner_sum_den * H_ch{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx);
                    
                    lambda_new{k_user}(1,idx)   = lambda_old{k_user}(1,idx) + (num/den);
                    lambda_old{k_user}(1,idx)   = lambda_new{k_user}(1,idx);
                    
                    V_new{k_user}(:,:,idx)  = pinv(sum_lambda{k_user}(:,:,idx) + lambda_new{k_user}(1,idx)*eye(Nt)) * H_ch{k_user, k_user}(:,:,idx)' * G_all{k_user}(:,:,idx);
                    V_norm{k_user}(:,:,idx) = norm(V_new{k_user}(:,:,idx), 2)^2;
                    
                    if V_norm{k_user}(:,:,idx) == 1
                        break;
                    end
                end
            end
        end
        
        for k_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = zeros(Nr, Nt);
            for j_user = 1:N_user
                sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_ch{k_user, j_user}(:,:,idx)*V_new{j_user}(:,:,idx)*V_new{j_user}(:,:,idx)'*H_ch{k_user, j_user}(:,:,idx)';
            end
            G{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + VarN(k_user)*eye(Nr)) * H_ch{k_user, k_user}(:,:,idx) * V_new{k_user}(:,:,idx);
            G_all{k_user}(:,:,idx) = G{k_user}(:,:,idx);
        end

        
%         for k_user = 1:N_user
%             for j_user = 1:N_user
%                 if k_user ~= j_user
%                     sum_MSE{k_user}(:,:,idx) = sum_MSE{k_user}(:,:,idx) + abs(G{k_user}(:,:,idx)' * H_ch{k_user, j_user}(:,:,idx) * V_new{j_user}(:,:,idx))^2;
%                 end
%             end
%             epslon{k_user}(j, idx) = abs((G{k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx) * V_new{k_user}(:,:,idx)) - 1)^2 + sum_MSE{k_user}(:,:,idx) + norm(G{k_user}(:,:,idx), 2)^2 * VarN(k_user);
%             eta_sum(j) = eta_sum(j)+ epslon{k_user}(j, idx); % Sum mean square error
%         end
        
        
        for k_user = 1:N_user
            Convergence_check(j, idx) = Convergence_check(j, idx) + norm(V_new{k_user}(:,:,idx) - V_old{k_user}(:,:,idx), 2);
        end
        if Convergence_check(j, idx) <= 10^-4
            break;
        end
    end
end
end