clear all
close all

I = double(pgmload('lenna.256.pgm'));

subplot(321)
imagesc(I),colormap(gray);

[nblig,nbcol] = size(I)
moy = mean(I(:))
I = I - moy;
It = bdct(I,[8 8]);

subplot(322)
imagesc(It);


	
ItDC = It(1,:);
It_hist = histc(ItDC,min(ItDC):10:max(ItDC));
ItDC_r = reshape(ItDC, nblig/8,nbcol/8);
subplot(323),imagesc(ItDC_r)
subplot(324),plot(min(ItDC):10:max(ItDC),It_hist)

% Variance de chaque sous-bande
It_var = var(It',1);
subplot(325),semilogy(It_var);


Qtables;

% Après balayage Zig-Zag
Itz = It(Zig,:);
Itz_var = var(Itz',1);
subplot(326),semilogy(Itz_var);
whos Itz

% Quantification
gamma = 0.5;
q = gamma*Quant(:); % Balayage des coefficients de quantification
for i=1:64
   It_quant(i,:) = fix(It(i,:)/q(i)); % arrondi vers 0 des coefficients
end

% Balayage Zig-Zag des coefficients quantifiés
It_quant_zig = It_quant(Zig,:);

% Compression
[bitstream,res] = JPEGlike(0, It_quant_zig);

% Reconstitution de l'image compressée
It_rec_zig = JPEGlike(bitstream, size(It_quant_zig,1), size(It_quant_zig,2));
It_rec(Zig,:) = It_rec_zig;

% Quantification inverse
for i = 1:64
    It_rec(i,:) = It_rec(i,:)*q(i);
end

I_rec = ibdct(It_rec, [8 8], [nblig nbcol]);

I_diff = I - I_rec;
SNR = 10*log10(var(I(:),1)/var(I_diff(:),1))
PSNR = 10*log10(255^2/var(I_diff(:),1))

figure(2)
subplot(231),colormap(gray)
imagesc(I);
subplot(232)
imagesc(I_rec);
subplot(233)
imagesc(I_diff);


