% Fonction de calcul de l'image pr�dite et du d�bit n�cessaire aux vecteurs mouvement
% vecteurs mouvement � 1 pixel avec une precision de +/-8 pixels 
%
% Entr�es:
%		compY : image originale (n)
%		compyYr : image reconstruite (n-1)
%		blksize : taille en pixels du bbloc sur lequel se fait la compensation en mouvement
%     SearchRange : port�e de la recherche des VM
%
% Sortie :
%		compPred : image pr�dite
%     MFx, MFy : champ de vecteurs de mouvement
%     bppMV : d�bit n�cessaire au codage des MV
%
% Michel Kieffer, le 31/01/2006

function [compPred,MFx,MFy,bppMV] = PicturePred(compY,compYr,blksize,SearchRange)

[nblig,nbcol] = size(compY);
bppMV = 0;

for lig = 1:blksize:nblig
   for col = 1:blksize:nbcol
      % Extraction des blocs de recherche dans l'image pr�c�dente
      SearchArea = 128*ones(blksize+2*SearchRange);
      if (lig==1) % Premi�re ligne de blocs
         if (col==1)
            SearchArea(SearchRange+1:2*SearchRange+blksize,...
               SearchRange+1:2*SearchRange+blksize) = compYr(lig:lig+blksize+SearchRange-1,...
               col:col+blksize+SearchRange-1);
         elseif (col==nbcol-blksize+1)
         	SearchArea(SearchRange+1:2*SearchRange+blksize,...
               1:SearchRange+blksize) = compYr(lig:lig+blksize+SearchRange-1,...
               col-SearchRange:col+blksize-1);
         else
         	SearchArea(SearchRange+1:2*SearchRange+blksize,...
               1:2*SearchRange+blksize) = compYr(lig:lig+blksize+SearchRange-1,...
               col-SearchRange:col+blksize+SearchRange-1);
         end
      elseif (lig==nblig-blksize+1) % derni�re ligne de blocs
         if (col==1)
            SearchArea( 1:SearchRange+blksize,...
               SearchRange+1:2*SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize-1,...
               col:col+blksize+SearchRange-1);
         elseif (col==nbcol-blksize+1)
         	SearchArea( 1:SearchRange+blksize,...
               1:SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize-1,...
               col-SearchRange:col+blksize-1);
         else
         	SearchArea( 1:SearchRange+blksize,...
               1:2*SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize-1,...
               col-SearchRange:col+blksize+SearchRange-1);
         end
      else % Autres lignes de blocs
         if (col==1)
            SearchArea(1:2*SearchRange+blksize,...
               SearchRange+1:2*SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize+SearchRange-1,...
               col:col+blksize+SearchRange-1);
         elseif (col==nbcol-blksize+1)
         	SearchArea(1:2*SearchRange+blksize,...
               1:SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize+SearchRange-1,...
               col-SearchRange:col+blksize-1);
         else
         	SearchArea(1:2*SearchRange+blksize,...
               1:2*SearchRange+blksize) = compYr(lig-SearchRange:lig+blksize+SearchRange-1,...
               col-SearchRange:col+blksize+SearchRange-1);
         end
      end
      % Recherche du vecteur mouvement
      [MVx,MVy,MCBloc] = MVSearch(compY(lig:lig+blksize-1,col:col+blksize-1),SearchArea,SearchRange);
      
      MFx(1+(lig-1)/blksize,1+(col-1)/blksize)=MVx;
      MFy(1+(lig-1)/blksize,1+(col-1)/blksize)=MVy;
      compPred(lig:lig+blksize-1,col:col+blksize-1) = MCBloc;
      
   end
end            
            
% Calcul du d�bit n�cessaire au codage des vecteurs mouvement
% Utilisation d'un codage similaire � celui des coefficients DC par JPEG
[bits_x,Res_x] = MVcode(1, MFx);
[bits_y,Res_y] = MVcode(1, MFy);
bppMV = (sum(Res_x([2:4])) + sum(Res_x([2:4])))/(nblig*nbcol);


