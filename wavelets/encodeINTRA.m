function [I_rec,PSNR,bpp] = encodeINTRA(I,gamma)

[nblig,nbcol] = size(I);
moy = mean(I(:));
I = I - moy;
It = bdct(I,[8 8]);

% Chargement de la table de quantification
Qtables;

% Après balayage Zig-Zag
Itz = It(Zig,:);

% Quantification
q = gamma*Quant(:); % Balayage des coefficients de quantification
for i=1:64
   It_quant(i,:) = fix(It(i,:)/q(i)); % arrondi vers 0 des coefficients
end

% Balayage Zig-Zag des coefficients quantifiés
It_quant_zig = It_quant(Zig,:);

% Compression
[bitstream,res] = JPEGlike(0, It_quant_zig);
bpp = length(bitstream)/(nblig*nbcol)

% Reconstitution de l'image compressée
It_rec_zig = JPEGlike(bitstream, size(It_quant_zig,1), size(It_quant_zig,2));
It_rec(Zig,:) = It_rec_zig;

% Quantification inverse
for i = 1:64
    It_rec(i,:) = It_rec(i,:)*q(i);
end

I_rec = ibdct(It_rec, [8 8], [nblig nbcol]) + moy;

I_diff = I - I_rec;
SNR = 10*log10(var(I(:),1)/var(I_diff(:),1));
PSNR = 10*log10(255^2/var(I_diff(:),1));
