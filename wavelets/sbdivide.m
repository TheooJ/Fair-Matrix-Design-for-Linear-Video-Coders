function Y=sbdivide(X, S)
%sbdivide - Split an array into separate subbands
%------------------------------------------------------------------------------
%SYNOPSIS       Y = sbdivide(X, S)
%                  Divide a matrix containing a subband split image into
%                  separate subbands.
%                  X is the image to divide.
%                  S is the number of splits.
%                  The result is a cell array with 3*S+1 elements, each
%                  containing one subband.
%
%                  The values in Y are ordered starting with the lowpass band.
%                  For instance, for a 2 level split, the values are in the
%                  following order:
%
%                            +------+------+-------------+
%                            |      |      |             |
%                            |  1   |  2   |             |
%                            |      |      |             |
%                            |------+------+      5      |
%                            |      |      |             |
%                            |  3   |  4   |             |
%                            |      |      |             |
%                            +------+------+-------------+
%                            |             |             |
%                            |             |             |
%                            |             |             |
%                            |      6      |      7      |
%                            |             |             |
%                            |             |             |
%                            |             |             |
%                            +-------------+-------------+
%
%SEE ALSO	sbmerge, dsbt2, idsbt2
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2003 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: sbdivide.m,v 1.2 2003/11/17 16:42:46 harna Exp harna $

dim=size(X);

if S==0
  Y={X};
else
  if min(dim)==1 %One-dimensional signal
    Y={sbdivide(X(1:ceil(length(X)/2)), S-1) {X(ceil(length(X)/2)+1:end)}};
  else
    r=dim(1);
    c=dim(2);
    Y=sbdivide(X(1:ceil(r/2),1:ceil(c/2)), S-1);
    Y{end+1} = X(1:ceil(r/2), ceil(c/2)+1:c);
    Y{end+1} = X(ceil(r/2)+1:r, 1:ceil(c/2));
    Y{end+1} = X(ceil(r/2)+1:r, ceil(c/2)+1:c);
  end
end
