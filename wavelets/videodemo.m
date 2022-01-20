clear all
close all

% Compression pour les images INTRA
gammaINTRA = 1;

% Taille des blocs pour la compensation du mouvement
blksize = 16;

% Précision de la recherche
searchRange = 8;

fid = fopen('foreman.qcif','r');
[compY,compU,compV]=yuv_readimage(fid,'qcif');

subplot(2,5,1)
imagesc(compY),colormap(gray)

%% PREMIERE PARTIE : Codage INTRA DE TOUTES LES TRAMES

% Codage INTRA de la première Image
[compYr,PSNR,nbbits] = encodeINTRA(compY,gammaINTRA);
figure(1),subplot(2,5,6)
imagesc(compYr),colormap(gray)
PSNR
nbbits

% Lecture et codage INTRA des images suivantes
for i=2:5
	[compY,compU,compV]=yuv_readimage(fid,'qcif');
	subplot(2,5,i)
	imagesc(compY),colormap(gray)
   [compYr,PSNR(i),nbbits(i)] = encodeINTRA(compY,gammaINTRA);
	subplot(2,5,5+i)
	imagesc(compYr),colormap(gray)
   [PSNR;nbbits]   
end

pause

%% SECONDE PARTIE : Codage INTRA DE LA PREMIERE TRAME ET INTER DES SUIVANTES
%%                  PAS DE COMPENSATION EN MOUVEMENT

% Codage INTRA de la première Image
[compYr,PSNR,nbbits] = encodeINTRA(compY,gammaINTRA);
figure(2),subplot(2,5,6)
imagesc(compYr),colormap(gray)
PSNR
nbbits

% Lecture et codage INTRA des images suivantes
for i=2:5
   [compY,compU,compV]=yuv_readimage(fid,'qcif');
	subplot(2,5,i)
	imagesc(compY),colormap(gray)
   
   % Calcul de la différence
   compDiff = compY - compYr;
   
   % Codage de la différence
   [compDiffr,PSNR(i),nbbits(i)] = encodeINTRA(compDiff,gammaINTRA);
   
   % Calcul de la nouvelle image reconstruite
   compYr = compYr + compDiffr;
   
   subplot(2,5,5+i)
	imagesc(compYr),colormap(gray)
   [PSNR;nbbits]   
end

pause

%% TROISIEME PARTIE : Codage INTRA DE LA PREMIERE TRAME ET INTER DES SUIVANTES
%%                    COMPENSATION EN MOUVEMENT

% Codage INTRA de la première Image
[compYr,PSNR,nbbits] = encodeINTRA(compY,gammaINTRA);
figure(3),subplot(2,5,6)
imagesc(compYr),colormap(gray)
PSNR
nbbits

% Lecture et codage INTRA des images suivantes
for i=2:5
   [compY,compU,compV]=yuv_readimage(fid,'qcif');
	subplot(2,5,i)
	imagesc(compY),colormap(gray)
   
   [compPred,MVx(:,:,i),MVy(:,:,i),nbbitsMV(i)] = PicturePred(compY,compYr,blksize,searchRange);
   
   % Calcul de la différence
   compDiff = compY - compPred;
   
   % Codage de la différence
   [compDiffr,PSNR(i),nbbits(i)] = encodeINTRA(compDiff,gammaINTRA);
   
   % Calcul de la nouvelle image reconstruite
   compYr = compPred + compDiffr;
   
   subplot(2,5,5+i)
	imagesc(compYr),colormap(gray)
   [PSNR;nbbits;nbbitsMV]   
end
