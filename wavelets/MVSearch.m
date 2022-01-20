% Recherche du vecteur mouvement minimisant la SAD
%
% Entrées : 
%	BlocY		 	: bloc à compenser en mouvement
%  SearchArea	: espace de recherche
%  SearchRange : largeur de la recherche
%
% Sorties :
%  MVx,MVy		: vecteurs mouvement
%  MCBloc      : bloc compensé en mouvement
%
% Michel Kieffer, le 31/01/2006

function [MVx,MVy,MCBloc] = MVSearch(BlocY,SearchArea,SearchRange)

[blkdim_x,blkdim_y] = size(BlocY);
SAD_min = 10e10;

for lig=-SearchRange:SearchRange
   for col=-SearchRange:SearchRange
      rangel = SearchRange+1+lig:SearchRange+lig+blkdim_x;
      rangec = SearchRange+1+col:SearchRange+col+blkdim_y;
      
      SAD_cur = sum(sum(abs(BlocY-SearchArea(rangel,rangec))));
      if (SAD_cur<=SAD_min)
         SAD_min=SAD_cur;
         MVx = lig;
         MVy = col;
         MCBloc = SearchArea(rangel,rangec);
      end
   end
end

