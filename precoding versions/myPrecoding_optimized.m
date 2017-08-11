function [ V, G ] = myPrecoding_optimized(H_ch, VarN, N)
%PRECODING Summary of this function goes here
%   Detailed explanation goes here
% Date: 17/6/2017
%% Variables Initialization
N_user  = size(H_ch, 1);
Nr      = size(H_ch{1},1);
Nt      = size(H_ch{1},2);

G       = cell(N_user, 1);
V       = cell(N_user, 1);
v_init  = cell(N_user, 1);
V_init  = cell(N_user, 1);
V_norm  = cell(N_user, 1);

sum_G  = cell(N_user, 1);
sum_V  = cell(N_user, 1);
lambda = cell(N_user, 1);

count = 5;
%% Initializing v_MMSE_k & g_MMSE_k

for i_user = 1:N_user
    lambda{i_user} =  zeros(1, 1, N);
    v_init{i_user} =  randn(Nt, 1)+ 1i* randn(Nt, 1);
    V_init{i_user} =  fft(v_init{i_user}, N, 3)./sqrt(N);
    V{i_user}      =  zeros(Nt, 1, N);
    
    sum_G{i_user}  = zeros(Nr, Nt, N);
    
    V_init_cell    = squeeze(num2cell(V_init{i_user},[1 2]));
    V_init_cell_h  = squeeze(num2cell(conj(permute(V_init{i_user}, [2, 1, 3])),[1 2]));
    
    for j_user = 1:N_user
        H_cell_ij      = squeeze(num2cell(H_ch{i_user, j_user},[1 2]));
        H_cell_ij_h    = squeeze(num2cell(conj(permute(H_ch{i_user, j_user}, [2, 1, 3])),[1 2]));
        
        inner_sum_mul  = cellfun(@(a, b, c, d) a*b*c*d, H_cell_ij, V_init_cell, V_init_cell_h, H_cell_ij_h, 'UniformOutput', false);
        
        inner_sum_G    = cat(3, inner_sum_mul{:});
        sum_G{i_user}  = sum_G{i_user} + inner_sum_G;
    end
    
    inv_mat_G = zeros(Nr, Nt, N);
    for idx = 1:N
        inv_mat_G(:,:,idx) = pinv(sum_G{i_user}(:,:,idx) + VarN(i_user)*eye(Nr));
    end
    
    H_cell_ii    = squeeze(num2cell(H_ch{i_user, i_user},[1 2]));
    inv_mat_cell = squeeze(num2cell(inv_mat_G,[1 2]));
    
    G_mul        = cellfun(@(a, b, c) a*b*c, inv_mat_cell, H_cell_ii, V_init_cell, 'UniformOutput', false);
    G{i_user}    = cat(3, G_mul{:});
end


%% Update
for i_user = 1:N_user
    sum_V{i_user} = zeros(Nr, Nt, N);
    
    for j_user = 1:N_user
        H_cell_ij     = squeeze(num2cell(H_ch{i_user, j_user},[1 2]));
        G_cell        = squeeze(num2cell(G{i_user},[1 2]));
        G_cell_h      = squeeze(num2cell(conj(permute(G{i_user}, [2, 1, 3])),[1 2]));
        H_cell_ij_h   = squeeze(num2cell(conj(permute(H_ch{i_user, j_user}, [2, 1, 3])),[1 2]));
        inner_sum_mul = cellfun(@(a, b, c, d) a*b*c*d, H_cell_ij_h, G_cell, G_cell_h, H_cell_ij, 'UniformOutput', false);
        
        inner_sum_V   = cat(3, inner_sum_mul{:});
        sum_V{i_user} = sum_V{i_user} + inner_sum_V;
    end
    
    
    inv_mat_init = zeros(Nr, Nt, N);
    for idx = 1:N
        inv_mat_init(:,:,idx) = pinv(sum_V{i_user}(:,:,idx));
    end
    inv_mat_cell    = squeeze(num2cell(inv_mat_init,[1 2]));
    H_cell_ii       = squeeze(num2cell(H_ch{i_user, i_user},[1 2]));
    H_cell_ii_h     = squeeze(num2cell(conj(permute(H_ch{i_user, i_user}, [2, 1, 3])),[1 2]));
    V_init_mul      = cellfun(@(a, b, c) a*b*c, inv_mat_cell, H_cell_ii_h, G_cell, 'UniformOutput', false);
    V_init{i_user}  = cat(3, V_init_mul{:});
    
    V_norm_cell = cellfun(@(a) norm(a, 2)^2, V_init_cell, 'UniformOutput', false);
    V_norm{i_user} = cat(3, V_norm_cell{:});
    
    idx_norm = find(V_norm{i_user} <= 1);
    if ~isempty(idx_norm)
        temp                    = V_init{i_user}(:,:,idx_norm);
        V{i_user}(:,:,idx_norm) = temp;
    end
    
    %% Updating lambda
    
    inv_mat_V = zeros(Nr, Nt, N);
    for idx = 1:N
        inv_mat_V(:,:,idx) = pinv(sum_V{i_user}(:,:,idx));
    end
    
    for i = 1:count
        inv_mat_V_cell  = squeeze(num2cell(inv_mat_V,[1 2]));
        inv_mat_nom_mul = cellfun(@(a, b) a*b, inv_mat_V_cell, inv_mat_V_cell, 'UniformOutput', false);
        inv_mat_den_mul = cellfun(@(a, b, c) a*b*c, inv_mat_V_cell, inv_mat_V_cell, inv_mat_V_cell, 'UniformOutput', false);
                
        nom_mul = cellfun(@(a, b, c, d, e) a*b*c*d*e - 1, G_cell_h, H_cell_ii, inv_mat_nom_mul, H_cell_ii_h, G_cell, 'UniformOutput', false);
        den_mul = cellfun(@(a, b, c, d, e) 2*a*b*c*d*e, G_cell_h, H_cell_ii, inv_mat_den_mul, H_cell_ii_h, G_cell, 'UniformOutput', false);
        
        lambda_old       = squeeze(num2cell(lambda{i_user},[1 2]));
        
        lambda_new_cell  = cellfun(@(a, b, c) a + b/c, lambda_old, nom_mul, den_mul, 'UniformOutput', false);
        lambda{i_user}   = cat(3, lambda_new_cell{:});
        
        lambda{i_user}(lambda{i_user} < 0) = 0;
        
        %% Updating v_MMSE_k
        inv_mat_V = zeros(Nr, Nt, N);
        for idx = 1:N
            inv_mat_V(:,:,idx) = pinv(sum_V{i_user}(:,:,idx) + lambda{i_user}(:,:,idx)*eye(Nr));
        end
        inv_mat_cell_u = squeeze(num2cell(inv_mat_V,[1 2]));
        V_update_mul   = cellfun(@(a, b, c) a*b*c, inv_mat_cell_u, H_cell_ii_h, G_cell, 'UniformOutput', false);
        V{i_user}      = cat(3, V_update_mul{:});
        
        
        V_norm_cell    = cellfun(@(a) norm(a, 2)^2, V_init_cell, 'UniformOutput', false);
        V_norm{i_user} = cat(3, V_norm_cell{:});
        
        % doesn't make sense
        if ~isempty(V_norm{i_user}(V_norm{i_user} == 1))
            break;
        end
    end
    %% Updating g_MMSE_k
    
    sum_G{i_user} = zeros(Nr, Nt, N);
    V_cell        = squeeze(num2cell(V{i_user},[1 2]));
    V_cell_h      = squeeze(num2cell(conj(permute(V{i_user}, [2, 1, 3])),[1 2]));
    
    for j_user = 1:N_user
        H_cell_ij      = squeeze(num2cell(H_ch{i_user, j_user},[1 2]));
        H_cell_ij_h    = squeeze(num2cell(conj(permute(H_ch{i_user, j_user}, [2, 1, 3])),[1 2]));
        inner_sum_mul  = cellfun(@(a, b, c, d) a*b*c*d, H_cell_ij, V_cell, V_cell_h, H_cell_ij_h, 'UniformOutput', false);
        inner_sum_G    = cat(3, inner_sum_mul{:});
        
        sum_G{i_user} = sum_G{i_user} + inner_sum_G;
    end
    
    inv_mat_G = zeros(Nr, Nt, N);
    for idx = 1:N
        inv_mat_G(:,:,idx) = pinv(sum_G{i_user}(:,:,idx) + VarN(i_user)*eye(Nr));
    end
    
    inv_mat_cell = squeeze(num2cell(inv_mat_G,[1 2]));
    
    G_mul        = cellfun(@(a, b, c) a*b*c, inv_mat_cell, H_cell_ii, V_cell, 'UniformOutput', false);
    G{i_user}    = cat(3, G_mul{:});
end

end