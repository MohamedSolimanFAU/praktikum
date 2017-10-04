function [ rx_sc ] = scfdma_ch( tx_sc, h_ch, VarN, noise_type)
%SCFDMA_CH Summary of this function goes here
%   Detailed explanation goes here

N_user = size(h_ch, 2);
Nr     = size(h_ch{1}, 1);
qh     = size(h_ch{1}, 3);

rx_fc = zeros(Nr, length(tx_sc) + qh - 1);
for i_user = 1:N_user
    rx_fc  = rx_fc + convLin(tx_sc , h_ch{i_user});
end

% rx_fc  = convLin(tx_sc , h_ch);
switch noise_type
    case 'awgn'
        noise  = sqrt(0.5*VarN).*(randn(Nr,length(rx_fc))+ 1i*randn(Nr,length(rx_fc)));
    case 'noiseless'
        noise  = 0;
    otherwise
        noise  = sqrt(0.5*VarN).*(randn(Nr,length(rx_fc))+ 1i*randn(Nr,length(rx_fc)));
end
rx_sc  = rx_fc  + noise ;

end