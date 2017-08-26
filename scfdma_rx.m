function [ rx ] = scfdma_rx( rx_sc, Scfdma, G )
%SCFDMA_RX Summary of this function goes here
%   Detailed explanation goes here

l_cp      = Scfdma.l_cp;
M         = Scfdma.M;
N         = Scfdma.N;
nu_0      = Scfdma.nu_0;
N_scSymb  = Scfdma.N_scSymb;

N_x = Scfdma.N + Scfdma.l_cp;

Nr      = size(G, 1);
rx_x    = zeros(Nr, N*N_scSymb);
Rx_x    = zeros(Nr, N*N_scSymb);
Rx_fir  = zeros(1, N*N_scSymb);
Rx      = zeros(1, M*N_scSymb);
rx      = zeros(1, N_scSymb*M);

rx_cp   = rx_sc;

% G' in switch case will be calculated
% NG      = size(G, 1);
% switch NG
%     case 1
%         G = squeeze(G)';
%     case 2
%         G = squeeze(permute(conj(G),[2 1 3]));
%     otherwise
%         G = squeeze(permute(conj(G),[2 1 3]));
% end

for i_bl = 0:N_scSymb-1
    rx_x(:, i_bl*N+(1:N))   = rx_cp(:, i_bl*N_x + (l_cp + 1:N_x));

    Rx_x(:, i_bl*N + (1:N)) = fft(rx_x(:, i_bl*N + (1:N)), N, 2)./sqrt(N);
    
    for idx = 1:N
        Rx_fir(:, i_bl*N + idx) = G(:,:,idx)' * Rx_x(:, i_bl*N + idx);
    end
    
    Rx(:, i_bl*M + (1:M))   = Rx_fir(:, i_bl*N + nu_0 +(1:M));

    rx(:, i_bl*M + (1:M))   = ifft(Rx(:, i_bl*M + (1:M)), M, 2).* sqrt(M);
end
end