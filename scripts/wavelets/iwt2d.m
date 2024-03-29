function out = iwt2d(in,code,stages)

% Usage: out = iwt2d(in,code,stages)
%
% Example for Daubechies D4
% in= wt2d([ 1 2 3 4 ; 5 6 7 8 ; 9 10 11 12 ; 13 14 15 16],'daub4',1); 
% code='daub4'; stages=1;
%
% Copyright (C) 2000-2001 F. Arguello (http://web.usc.es/~elusive)
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

  [N2 N1]=size(in);
  [h g hh gg desp]=wcoeff(code);

% EXAMPLE FOR DAUBECHES D4
% hh= [-0.129410, 0.224144, 0.836516, 0.482963];
% gg= [-0.482963, 0.836516, -0.224144, -0.129410];

  out=in;
  vsize=N2/(2^stages); hsize=N1/(2^stages);
  for k=1:stages
    hsize=hsize*2; vsize=vsize*2;
    out([1:vsize],[1:hsize]) ...
        = wavelet_inv_step(out([1:vsize],[1:hsize]),hh,gg,desp);
    out([1:vsize],[1:hsize]) ...
        = wavelet_inv_step(out([1:vsize],[1:hsize])',hh,gg,desp)';
  end
  
end
