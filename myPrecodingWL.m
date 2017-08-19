function [ V, G ] = myPrecodingWL( H_all, H_ch, VarN, N, WL )
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

Dr    = WL.Dr;
Br    = WL.Br;
Rs    = cell(N_user, 1);
%% Initializing Rs in frequency domain

for k_user = 1:N_user
    Rs{k_user} = fft(WL.Rs(:,:,k_user), N, 3);
end

%% Initializing beamforming vector
for idx = 1:N
    for k_user = 1:N_user
        [vector, value] = eig(H_ch{k_user, k_user}(:,:, idx)' * H_ch{k_user, k_user}(:,:, idx));
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
%% Initializing v_MMSE_k and g_MMSE_k

for k_user = 1:N_user
%     v_init{k_user} = randn(Nt, 1)+ 1i* randn(Nt, 1);
%     V_init{k_user} = fft(v_init{k_user}, N, 3);
    V{k_user}      = zeros(Dr, Br, N);
    G{k_user}      = zeros(2*Nt, Br, N);
    G_init{k_user} = zeros(2*Nt, Br, N);
    sum_G{k_user}  = zeros(2*Nr, 2*Nt, N);
    sum_V{k_user}  = zeros(2*Nr, 2*Nt, N);
end

%% update v_MMSE_k and g_MMSE_k
for idx = 1:N
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_all{k_user, j_user}(:,:,idx)* V_init{j_user}(:,:,idx)* Rs{j_user}(:,:,idx) * V_init{j_user}(:,:,idx)' * H_all{k_user, j_user}(:,:,idx)';
        end
        G_init{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + (VarN(k_user)/2)*eye(2*Nr)) * H_all{k_user, k_user}(:,:,idx) * V_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
    end
    
    for k_user = 1:N_user
        for j_user = 1:N_user
            sum_V{k_user}(:,:,idx)  = sum_V{k_user}(:,:,idx) + H_all{j_user, k_user}(:,:,idx)' * G_init{j_user}(:,:,idx) * Rs{j_user}(:,:,idx) * G_init{j_user}(:,:,idx)' * H_all{j_user, k_user}(:,:,idx);
        end
        V_init{k_user}(:,:,idx) = pinv(sum_V{k_user}(:,:,idx)) * H_all{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
        V_norm{k_user}(:,:,idx) = trace(V_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx) * V_init{k_user}(:,:,idx)');
    end
    % Kolo tamam l7d hena
    % start modify the next bit of code
    for k_user = 1:N_user
        if V_norm{k_user}(:,:,idx) <= 1
            V{k_user}(:,:,idx) = V_init{k_user}(:,:,idx);
        elseif V_norm{k_user}(:,:,idx) > 1
            lambda{k_user}(1,idx) =  0;
            
            for i = 1:count
                inverse_mat   = pinv(sum_V{k_user}(:,:,idx) + lambda{k_user}(1,idx)*eye(2*Nt));
                inner_sum_num = inverse_mat * inverse_mat;
                inner_sum_den = inverse_mat * inverse_mat * inverse_mat;
                
                num = trace( Rs{k_user}(:,:,idx) * G_init{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx)* inner_sum_num * H_all{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx)) - 1;
                den = trace( 2 * Rs{k_user}(:,:,idx) * G_init{k_user}(:,:,idx)' * H_all{k_user, k_user}(:,:,idx) * inner_sum_den * H_all{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx));
                
                lambda{k_user}(1,idx)   = lambda{k_user}(1,idx) + (num/den);
                
                V{k_user}(:,:,idx)      = pinv(sum_V{k_user}(:,:,idx) + lambda{k_user}(1,idx)*eye(2*Nt)) * H_all{k_user, k_user}(:,:,idx)' * G_init{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
                V_norm{k_user}(:,:,idx) = trace(V{k_user}(:,:,idx) *  Rs{k_user}(:,:,idx) * V{k_user}(:,:,idx)');
                
                if round(V_norm{k_user}(:,:,idx), 4) == 1
                    break;
                end
            end
        end
    end
    
    for k_user = 1:N_user
        sum_G{k_user}(:,:,idx)  = zeros(2*Nr, 2*Nt);
        for j_user = 1:N_user
            sum_G{k_user}(:,:,idx)  = sum_G{k_user}(:,:,idx) + H_all{k_user, j_user}(:,:,idx) * V{j_user}(:,:,idx) * Rs{j_user}(:,:,idx) * V{j_user}(:,:,idx)'*H_all{k_user, j_user}(:,:,idx)';
        end
        G{k_user}(:,:,idx) = pinv(sum_G{k_user}(:,:,idx) + (VarN(k_user)/2)*eye(2*Nr)) * H_all{k_user, k_user}(:,:,idx) * V{k_user}(:,:,idx) * Rs{k_user}(:,:,idx);
    end
end

end