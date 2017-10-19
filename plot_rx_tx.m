uuser =3;

figure;
subplot(2,1,1)
stem(real(rx{uuser}),'rx')
hold on
stem(real(tx{uuser}),'bo')

subplot(2,1,2)
stem(imag(rx{uuser}),'rx')
hold on
stem(imag(tx{uuser}),'bo')
