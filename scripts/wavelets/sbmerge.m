function X=sbmerge(Y)
%sbmerge - Merge a cellarray of subbands into a single array
%------------------------------------------------------------------------------
%SYNOPSIS       Y = sbmerge(X)
%                  Merge a cell array of subbands into a single matrix.
%
%                  The elements of X are assumed to be ordered starting with
%                  the lowpass band.
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
%SEE ALSO	sbdivide, dsbt2, idsbt2
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2003 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: sbmerge.m,v 1.1 2003/11/17 16:23:31 harna Exp $

dim=size(Y{end});

if min(dim)==1 %One-dimensional signal
  if length(Y)==1
    X=Y{1};
  else
    X=[sbmerge(Y(1:end-1)) Y{end}];
  end
else
  if length(Y)==1
    X=Y{1};
  else
    X=[sbmerge(Y(1:end-3)) Y{end-2}; Y{end-1} Y{end}];
  end
end
