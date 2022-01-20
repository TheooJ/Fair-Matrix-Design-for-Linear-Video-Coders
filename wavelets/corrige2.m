clear all
close all

I = double(pgmload('lenna.256.pgm'));

I = (I-mean(I(:)))/256;

Ito = wt2d(I,'filter9-7',2);

subplot(221)
imagesc(I),colormap(gray);

subplot(222)
imagesc(Ito),colormap(gray);

Itod = sbdivide(Ito,2); % 2 correspond au nombre de niveaux
for i=1:length(Itod)
%    figure, imagesc(Itod{i}),colormap(gray)
end

Q2 = [0.1 0.2 0.2 0.2 0.4 0.4 0.4];
Itodq2 = sbquant(Itod,Q2);
Itodr2 = sbrec(Itodq2,Q2);

% Regroupement des sous-bandes
Itom = sbmerge(Itodr2);

I_hat = iwt2d(Itom,'filter9-7',2);

10*log10(mean((I(:)).^2)/mean((I(:)-I_hat(:)).^2))

subplot(223)
imagesc(I_hat),colormap(gray);
