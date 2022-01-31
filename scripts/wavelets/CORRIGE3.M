clear all
close all

I = double(pgmload('lenna.256.pgm'));

%I = (I-mean(I(:)))/256;

Ito = wt2d(I,'filter9-7',1);

Itod = sbdivide(Ito,1); % 2 correspond au nombre de niveaux
Itod{2} = 3*Itod{2}+256;
Itod{3} = 3*Itod{3}+256;
Itod{4} = 3*Itod{4}+256;
% Itod{5} = 3*Itod{5}+1024;
% Itod{6} = 3*Itod{6}+1024;
% Itod{7} = 3*Itod{7}+1024;
% Itod{8} = 3*Itod{8}+1024;
% Itod{9} = 3*Itod{9}+1024;
% Itod{10} = 3*Itod{10}+1024;
Itom = sbmerge(Itod);


subplot(221)
imagesc(I),colormap(gray);

subplot(222)
imagesc(Itom,[0,511]),colormap(gray);

for i=1:length(Itod)
%    figure, imagesc(Itod{i}),colormap(gray)
end

Q2 = [0.1 0.2 0.2 0.2];
Itodq2 = sbquant(Itod,Q2);
Itodr2 = sbrec(Itodq2,Q2);

% Regroupement des sous-bandes
Itom = sbmerge(Itodr2);

I_hat = iwt2d(Itom,'filter9-7',3);

10*log10(mean((I(:)).^2)/mean((I(:)-I_hat(:)).^2))

subplot(223)
imagesc(I_hat),colormap(gray);
