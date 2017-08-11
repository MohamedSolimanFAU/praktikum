function [ v, g ] = myprecoding_3( H_ch, VarN, N )
%MYPRECODING_2 Summary of this function goes here
%   Detailed explanation goes here

%% Variables Initialization
N_user  = size(H_ch, 1);
Nr      = size(H_ch{1},1);
Nt      = size(H_ch{1},2);

g       = cell(N_user, 1);
v       = cell(N_user, 1);
v_init  = cell(N_user, 1);
v_norm  = cell(N_user, 1);

sum_g  = cell(N_user, 1);
sum_v  = cell(N_user, 1);
sum_l  = cell(N_user, 1);
lambda = cell(N_user, 1);

count = 5;

%% Initializing v_MMSE_k & g_MMSE_k

for i_user = 1:N_user
	v_init{i_user} = randn(Nt, 1, N)+ 1i* randn(Nt, 1, N);
    v{i_user}      = zeros(Nt, 1, N);
    g{i_user}      = zeros(Nt, 1, N);
    sum_g{i_user}  = zeros(Nr, Nt, N);
    sum_v{i_user}  = zeros(Nr, Nt, N);
    sum_l{i_user}  = zeros(Nr, Nt, N);
end

for idx = 1:N
    for i_user = 1:N_user
        for j_user = 1:N_user
            sum_g{i_user}(:,:,idx)  = sum_g{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)*v_init{i_user}(:,:,idx)*v_init{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx)';
        end
        g{i_user}(:,:,idx) = pinv(sum_g{i_user}(:,:,idx) + VarN(i_user)*eye(Nr)) * H_ch{i_user, i_user}(:,:,idx) * v_init{i_user}(:,:,idx);
    end
    
    for i_user = 1:N_user
        for j_user = 1:N_user
            sum_v{i_user}(:,:,idx)  = sum_v{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)'*g{i_user}(:,:,idx)*g{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx);
        end
        v_init{i_user}(:,:,idx) = pinv(sum_v{i_user}(:,:,idx)) * H_ch{i_user, i_user}(:,:,idx)' * g{i_user}(:,:,idx);
        v_norm{i_user}(:,:,idx) = norm(v_init{i_user}(:,:,idx), 2)^2;
    end
    
    for i_user = 1:N_user
        if v_norm{i_user}(:,:,idx) <= 1
            v{i_user}(:,:,idx) = v_init{i_user}(:,:,idx);
        elseif v_norm{i_user}(:,:,idx) > 1
            for j_user = 1:N_user
                sum_l{i_user}(:,:,idx)  = sum_l{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)'*g{i_user}(:,:,idx)*g{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx);
            end
            
            lambda{i_user}(1,idx) =  0;
            
            for i = 1:count
                inverse_mat   = pinv(sum_l{i_user}(:,:,idx) + lambda{i_user}(1,idx)*eye(Nt));
                inner_sum_num = inverse_mat * inverse_mat;
                inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                
                num =     g{i_user}(:,:,idx)' * H_ch{i_user, i_user}(:,:,idx)* inner_sum_num * H_ch{i_user, i_user}(:,:,idx)' * g{i_user}(:,:,idx) - 1;
                den = 2 * g{i_user}(:,:,idx)' * H_ch{i_user, i_user}(:,:,idx)* inner_sum_den * H_ch{i_user, i_user}(:,:,idx)' * g{i_user}(:,:,idx);
                
                lambda{i_user}(1,idx) = lambda{i_user}(1,idx) + (num/den);
                
                v{i_user}(:,:,idx)      = pinv(sum_l{i_user}(:,:,idx) + lambda{i_user}(1,idx)*eye(Nt)) * H_ch{i_user, i_user}(:,:,idx)' * g{i_user}(:,:,idx);
                v_norm{i_user}(:,:,idx) = norm(v{i_user}(:,:,idx), 2)^2;
                
                if v_norm{i_user}(:,:,idx) == 1
                    break;
                end
            end
        end
    end
    
    for i_user = 1:N_user
        for j_user = 1:N_user
            sum_g{i_user}(:,:,idx)  = sum_g{i_user}(:,:,idx) + H_ch{i_user, j_user}(:,:,idx)*v{i_user}(:,:,idx)*v{i_user}(:,:,idx)'*H_ch{i_user, j_user}(:,:,idx)';
        end
        g{i_user}(:,:,idx) = pinv(sum_g{i_user}(:,:,idx) + VarN(i_user)*eye(Nr)) * H_ch{i_user, i_user}(:,:,idx) * v{i_user}(:,:,idx);
    end
    
    
end

end

