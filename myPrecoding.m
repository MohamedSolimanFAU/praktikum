function [ V, G ] = myPrecoding( H_ch, VarN, N )
%MYPRECODING_2 Summary of this function goes here
%   Detailed explanation goes here

%% Variables Initialization
N_user  = size(H_ch, 1);
Nr      = size(H_ch{1},1);
Nt      = size(H_ch{1},2);

G       = cell(N_user, 1);
G_init	= cell(N_user, 1);
V       = cell(N_user, 1);
v_init  = cell(N_user, 1);
V_init  = cell(N_user, 1);
V_norm  = cell(N_user, 1);

sum_G  = cell(N_user, 1);
sum_V  = cell(N_user, 1);
lambda = cell(N_user, 1);

count = 10;

%% Initializing v_MMSE_k & g_MMSE_k

for k_user = 1:N_user
    v_init{k_user} = randn(Nt, 1)+ 1i* randn(Nt, 1);
    V_init{k_user} = fft(v_init{k_user}, N, 3);
    V{k_user}      = zeros(Nt, 1, N);
    G{k_user}      = zeros(Nt, 1, N);
    G_init{k_user} = zeros(Nt, 1, N);
    sum_G{k_user}  = zeros(Nr, Nt, N);
    sum_V{k_user}  = zeros(Nr, Nt, N);
end

for idx = 1:N
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_ch{k_user, j_user}(:,:,idx)*V_init{j_user}(:,:,idx)*V_init{j_user}(:,:,idx)'*H_ch{k_user, j_user}(:,:,idx)';
        end
        G_init{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + VarN(k_user)*eye(Nr)) * H_ch{k_user, k_user}(:,:,idx) * V_init{k_user}(:,:,idx);
    end
    
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_V{k_user}(:,:,idx)  = sum_V{k_user}(:,:,idx) + H_ch{j_user, k_user}(:,:,idx)'*G_init{j_user}(:,:,idx)*G_init{j_user}(:,:,idx)'*H_ch{j_user, k_user}(:,:,idx);
        end
        V_init{k_user}(:,:,idx) = pinv(sum_V{k_user}(:,:,idx)) * H_ch{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx);
        V_norm{k_user}(:,:,idx) = norm(V_init{k_user}(:,:,idx), 2)^2;
    end
    
    for k_user = 1:N_user
        if V_norm{k_user}(:,:,idx) <= 1
            V{k_user}(:,:,idx) = V_init{k_user}(:,:,idx);
        elseif V_norm{k_user}(:,:,idx) > 1
            lambda{k_user}(1,idx) =  0;
            
            for i = 1:count
                inverse_mat   = pinv(sum_V{k_user}(:,:,idx) + lambda{k_user}(1,idx)*eye(Nt));
                inner_sum_num = inverse_mat * inverse_mat;
                inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                
                num =     G_init{k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx)* inner_sum_num * H_ch{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx) - 1;
                den = 2 * G_init{k_user}(:,:,idx)' * H_ch{k_user, k_user}(:,:,idx)* inner_sum_den * H_ch{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx);
                
                lambda{k_user}(1,idx)   = lambda{k_user}(1,idx) + (num/den);
                
                V{k_user}(:,:,idx)      = pinv(sum_V{k_user}(:,:,idx) + lambda{k_user}(1,idx)*eye(Nt)) * H_ch{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx);
                V_norm{k_user}(:,:,idx) = norm(V{k_user}(:,:,idx), 2)^2;
                
                if round(V_norm{k_user}(:,:,idx), 4) == 1
                    break;
                end
            end
        end
    end
    
    for k_user = 1:N_user
        sum_G{k_user}(:,:,idx)  = zeros(Nr, Nt);
        for j_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_ch{k_user, j_user}(:,:,idx)*V{j_user}(:,:,idx)*V{j_user}(:,:,idx)'*H_ch{k_user, j_user}(:,:,idx)';
        end
        G{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + VarN(k_user)*eye(Nr)) * H_ch{k_user, k_user}(:,:,idx) * V{k_user}(:,:,idx);
    end
end

end