clear all
close all

I = double(pgmload('lenna.256.pgm'));

I = I-mean(I(:));

subplot(221)
imagesc(I),colormap(gray);

It = bdct(I,[8 8]);
subplot(222)
imagesc(It);

Ite = It(2,:);

a = hist(Ite,min(Ite):max(Ite));

Ite_r = reshape(Ite,256/8,256/8);
subplot(223)
imagesc(Ite_r);
subplot(224)
plot(min(Ite):max(Ite),a);

It_var = var(It',1);
semilogy(It_var)


Qtables; % Charge le vecteur Zig qui permet de faire le balayage zig-zag
gamma = 4;
	
q = gamma*Quant(Zig);
	
for i=1:64
  It_quant(i,:) = fix(It(i,:)/q(i));
end

% Codage différentiel
It_quant(1,2:end) = It_quant(1,2:end)-It_quant(1,1:end-1);

for i=1:64
    ent(i) = entropy(It_quant(i,:));
end
entropy(It_quant(:))
mean(ent)
ent

for i=2:length(It_quant(1,:))
    It_quant(1,i) = It_quant(1,i)+It_quant(1,i-1);
end

for i=1:64
  It_rec(i,:) = It_quant(i,:)*q(i);
end

I_rec = ibdct(It_rec,[8 8],[256 256]);

subplot(224)
imagesc(I_rec);