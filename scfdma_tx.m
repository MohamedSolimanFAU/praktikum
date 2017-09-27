function [ tx_cp ] = scfdma_tx( tx, Scfdma, V )
%SCFDMA_TX Summary of this function goes here
%   Detailed explanation goes here

% V : [Nt x 1 x N]

l_cp      = Scfdma.l_cp;
M         = Scfdma.M;
N         = Scfdma.N;
nu_0      = Scfdma.nu_0;
N_scSymb  = Scfdma.N_scSymb;

N_x = Scfdma.N + Scfdma.l_cp;

Nt      = size(V, 1);
Tx      = zeros(1, M*N_scSymb);
Tx_x    = zeros(Nt, N*N_scSymb);
tx_x    = zeros(Nt, N*N_scSymb);
tx_cp   = zeros(Nt, N*N_scSymb);


for i_bl = 0:N_scSymb-1
    Tx(:, i_bl*M+(1:M))          = fft(tx(:, i_bl*M+(1:M)),M, 2)./sqrt(M);
    
    for idx = 1:M
        Tx_x(:, i_bl*N+nu_0+idx)   = V(:,:,idx+nu_0) * Tx(:, i_bl*M+idx);
    end
    
    tx_x(:, i_bl*N+(1:N))        = ifft(Tx_x(:, i_bl*N+(1:N)), N, 2).*sqrt(N);
    
    tx_cp(:, i_bl*N_x+(1:N_x))   = [tx_x(:, i_bl*N+(N-l_cp+1:N)) tx_x(:, i_bl*N+(1:N))];
end
end