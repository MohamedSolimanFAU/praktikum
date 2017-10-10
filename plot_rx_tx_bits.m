uuser =1;
le = 1:80;
figure;
title('bits')
stem(real(rx_bits{uuser}(le)),'rx')
hold on
stem(real(tx_bits{uuser}(le)),'bo')