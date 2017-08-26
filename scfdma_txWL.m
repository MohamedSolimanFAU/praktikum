function [ tx_cp ] = scfdma_txWL( tx, Scfdma, V )
%SCFDMA_TX Summary of this function goes here
%   Detailed explanation goes here

% V : [Dr x Br x N]

l_cp      = Scfdma.l_cp;
M         = Scfdma.M;
N         = Scfdma.N;
nu_0      = Scfdma.nu_0;
N_scSymb  = Scfdma.N_scSymb;

N_x = Scfdma.N + Scfdma.l_cp;

Nt      = size(V, 1);
Tx      = zeros(1, M*N_scSymb);
Tx_aug  = zeros(Nt/2, M*N_scSymb);
Tx_x    = zeros(Nt, N*N_scSymb);
Tx_com  = zeros(Nt/2, N*N_scSymb);
tx_x    = zeros(Nt/2, N*N_scSymb);
tx_cp   = zeros(Nt/2, N*N_scSymb);

V       = squeeze(V(:,:, nu_0+(1:M)));

for i_bl = 0:N_scSymb-1
    Tx(:, i_bl*M+(1:M))          = fft(tx(:, i_bl*M+(1:M)),M, 2)./sqrt(M);
    
    Tx_aug(:, i_bl*M+(1:M))      = [real(Tx(:, i_bl*M+(1:M))); imag(Tx(:, i_bl*M+(1:M)))];
    
    for idx = 1:M
        Tx_x(:, i_bl*N+nu_0+idx)   = V(:,:,idx) * Tx_aug(:, i_bl*M+idx);
    end
    % combine
    for idx = 1:N
        Tx_com(1, i_bl*N+idx) = Tx_x(1, i_bl*N+idx) + 1i*Tx_x(3, i_bl*N+idx);
        Tx_com(2, i_bl*N+idx) = Tx_x(2, i_bl*N+idx) + 1i*Tx_x(4, i_bl*N+idx);
    end
    
    tx_x(:, i_bl*N+(1:N))        = ifft(Tx_com(:, i_bl*N+(1:N)), N, 2).*sqrt(N);
    
    tx_cp(:, i_bl*N_x+(1:N_x))   = [tx_x(:, i_bl*N+(N-l_cp+1:N)) tx_x(:, i_bl*N+(1:N))];
end
end