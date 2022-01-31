% function [compY,compU,compV]=yuv_readimage(fid,format)
%
% Fonction de lexture d'une image (YUV) dans un fichier YUV
% Entree :
%    fid : identifiant d'un fichier (obtenu par fopen)
%    format ('cif' ou 'qcif')
%
% Sortie :
%    compY, compU et compV : composantes YUV de l'image
%
% Exemple :
%
% fid = fopen('foreman.qcif','r');
% [compY,compU,compV]=yuv_readimage(fid,'qcif')

function [compY,compU,compV]=yuv_readimage(fid,format)

% Calcul du format
if (strcmp(format,'cif'))
    width = 352;
    height = 288;
elseif (strcmp(format,'qcif'))
    width = 176;
    height = 144;
else
    width = 176;
    height = 144;
end
 
compY = fread(fid,[width height], 'uint8');
compY = compY';
compU = fread(fid, [width/2 height/2], 'uint8');
compV = fread(fid, [width/2 height/2], 'uint8');
 