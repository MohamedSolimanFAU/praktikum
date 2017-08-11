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
Tx      = zeros(Nt, M*N_scSymb);
Tx_x    = zeros(Nt, N*N_scSymb);
tx_x    = zeros(Nt, N*N_scSymb);
tx_cp   = zeros(Nt, N*N_scSymb);

tx_temp = tx;
V       = squeeze(V(:, nu_0+(1:M)));


for i_bl = 0:N_scSymb-1
    Tx(:, i_bl*M+(1:M))          = fft(tx_temp(:, i_bl*M+(1:M)),M, 2)./sqrt(M);
    
    Tx_x(:, i_bl*N+nu_0+(1:M))   = V .* Tx(:, i_bl*M+(1:M));
    
    tx_x(:, i_bl*N+(1:N))        = ifft(Tx_x(:, i_bl*N+(1:N)), N, 2).*sqrt(N);
    
    tx_cp(:, i_bl*N_x+(1:N_x))   = [tx_x(:, i_bl*N+(N-l_cp+1:N)) tx_x(:, i_bl*N+(1:N))];
end
end