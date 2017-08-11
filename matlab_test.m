%%
% How to transform a 3D matrix into a cell array, which will make the multiplication
% much easier.

% Mycell = num2cell(Mymatrix,[1 2]);
H_ch_cell      = squeeze(num2cell(H_ch{i_user, j_user},[1 2]));
V_init_cell    = squeeze(num2cell(V_init{i_user},[1 2]));
V_init_cell_t  = squeeze(num2cell(conj(permute(V_init{i_user}, [2, 1, 3])),[1 2]));
H_ch_cell_t    = squeeze(num2cell(conj(permute(H_ch{i_user, j_user}, [2, 1, 3])),[1 2]));

% To perform multiplication on this array use this command
C = cellfun(@(a, b, c, d) a*b*c*d, H_ch_cell, V_init_cell, V_init_cell_t, H_ch_cell_t, 'UniformOutput', false);

inv_mat = cellfun(@(a) pinv(a + VarN(i_user)*eye(Nr)), sum_G{i_user}, 'UniformOutput', false);

norm = cellfun(@(a) norm(a, 2)^2, V_init{i_user}, 'UniformOutput', false);
% Convert cell array back to 3D matrix
dim = 3;
cat(dim, C{:});


% Z = cellfun(@(x) reshape(x,Nr,Nt,[]), C, 'UniformOutput', 0);
% out = cell2mat(Z);
% user permute if the dimension is not relevant
