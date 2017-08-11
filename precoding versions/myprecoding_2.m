function [ V, G ] = myprecoding_2( H_ch, VarN, N )
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

for i_user = 1:N_user
    v_init{i_user} = randn(Nt, 1)+ 1i* randn(Nt, 1);
    V_init{i_user} = fft(v_init{i_user}, N, 3);
    V{i_user}      = zeros(Nt, 1, N);
    G{i_user}      = zeros(Nt, 1, N);
    G_init{i_user} = zeros(Nt, 1, N);
    sum_G{i_user}  = zeros(Nr, Nt, N);
    sum_V{i_user}  = zeros(Nr, Nt, N);
end

for idx = 1:N
    for i_user = 1:N_user
        for j_user = 1:N_user
            sum_G{i_user}(:,:,idx)  = sum_G{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)*V_init{i_user}(:,:,idx)*V_init{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx)';
        end
        G_init{i_user}(:,:,idx) = pinv(sum_G{i_user}(:,:,idx) + VarN(i_user)*eye(Nr)) * H_ch{i_user, i_user}(:,:,idx) * V_init{i_user}(:,:,idx);
    end
    
    for i_user = 1:N_user
        for j_user = 1:N_user
            sum_V{i_user}(:,:,idx)  = sum_V{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)'*G_init{i_user}(:,:,idx)*G_init{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx);
        end
        V_init{i_user}(:,:,idx) = pinv(sum_V{i_user}(:,:,idx)) * H_ch{i_user, i_user}(:,:,idx)' * G_init{i_user}(:,:,idx);
        V_norm{i_user}(:,:,idx) = norm(V_init{i_user}(:,:,idx), 2)^2;
    end
    
    for i_user = 1:N_user
        if V_norm{i_user}(:,:,idx) <= 1
            V{i_user}(:,:,idx) = V_init{i_user}(:,:,idx);
        elseif V_norm{i_user}(:,:,idx) > 1
            lambda{i_user}(1,idx) =  0;
            
            for i = 1:count
                inverse_mat   = pinv(sum_V{i_user}(:,:,idx) + lambda{i_user}(1,idx)*eye(Nt));
                inner_sum_num = inverse_mat * inverse_mat;
                inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                
                num =     G_init{i_user}(:,:,idx)' * H_ch{i_user, i_user}(:,:,idx)* inner_sum_num * H_ch{i_user, i_user}(:,:,idx)' * G_init{i_user}(:,:,idx) - 1;
                den = 2 * G_init{i_user}(:,:,idx)' * H_ch{i_user, i_user}(:,:,idx)* inner_sum_den * H_ch{i_user, i_user}(:,:,idx)' * G_init{i_user}(:,:,idx);
                
                lambda{i_user}(1,idx) = lambda{i_user}(1,idx) + (num/den);
                
                V{i_user}(:,:,idx)      = pinv(sum_V{i_user}(:,:,idx) + lambda{i_user}(1,idx)*eye(Nt)) * H_ch{i_user, i_user}(:,:,idx)' * G_init{i_user}(:,:,idx);
                V_norm{i_user}(:,:,idx) = norm(V{i_user}(:,:,idx), 2)^2;
                
                if round(V_norm{i_user}(:,:,idx), 4) == 1
                    break;
                end
            end
        end
    end
    
    for i_user = 1:N_user
        sum_G{i_user}(:,:,idx)  = zeros(Nr, Nt);
        for j_user = 1:N_user
            sum_G{i_user}(:,:,idx)  = sum_G{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)*V{i_user}(:,:,idx)*V{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx)';
        end
        G{i_user}(:,:,idx) = pinv(sum_G{i_user}(:,:,idx) + VarN(i_user)*eye(Nr)) * H_ch{i_user, i_user}(:,:,idx) * V{i_user}(:,:,idx);
    end
end

end

