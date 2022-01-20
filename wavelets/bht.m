function t=bht(im, M)
%bht - block hadamard transform of image
%------------------------------------------------------------------------------
%SYNOPSIS	t = bht(im, [M N])
%		  Perform blockwise hadamard trasnform on image im,
%                 using blocks of size MxN. The resulting transformed
%                 image t consists of one column per block.
%
%       	t = bht(im, M)
%                 As above, but with blocks of size MxM
%
%
%SEE ALSO	ibht, bdct, ibdct
%
%RCSID          $Id: bht.m,v 1.1 1998/11/22 11:37:49 harna Exp $
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 1998 Image Coding Group. LiU, SWEDEN

if (nargin == 0)
  error('No input arguments.')
end

if (nargin == 1)
  error('No blocksize given.')
end

if (length(M) == 1)
  M = [M M];
end

t=had2basemx(M)*imtocol(im, M, 'distinct');

