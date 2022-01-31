function tq=sbquant(t, q)
%sbquant - Do quantization of a subband transformed image
%------------------------------------------------------------------------------
%SYNOPSIS	tq = sbquant(t, q)
%		   where
%
%		     t    : image to quantize, either a single matrix or a
%                           cell array of subbands (such as produced by
%                           sbdivide).
%		     q    : Either a scalar, or a vector of quantization
%                           values, one for each subband of t. The length
%                           of q must be 3n+1, for a positive integer n.
%                           Note that you have top keep track of the number
%                           of subbands splits performed on t yourself, since
%                           there is no way for this function to know that.
%                           The values of q are step sizes for uniform
%                           quantizers.
%                           If q is a scalar, the same quantization step
%                           is used in all subbands.
%
%                           The values in q are ordered starting with the
%                           quantization value for the lowpass band.
%                           For instance, for a 2 level split, the values
%                           are in the following order:
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
%
%SEE ALSO	sbrec, dsbt2, idsbt2, sbdivide, sbmerge
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2001 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: sbquant.m,v 1.3 2003/11/17 16:42:04 harna Exp $

  
if (nargin ~= 2)
  error('Wrong number of input arguments.')
end

q = q(:);  % Just for convenience

if iscell(t)
  tq=t;
  if length(q)==1
    q=q*ones(length(tq));
  elseif length(q)~=length(t)
    error('Wrong number of quantization values.')
  end
  for k=1:length(t)
    tq{k}=round(tq{k}/q(k));
  end
else
  nq = floor(length(q)/3);

  if mod(length(q), 3) ~= 1
    error('Wrong number of quantization values.')
  end

  if length(q) == 1
    tq=round(t/q);   % Only one quantization value
  else
    tq=zeros(size(t));
    sy=floor(size(t, 1)/(2^nq));
    sx=floor(size(t, 2)/(2^nq));
    tq(1:sy,1:sx) = round(t(1:sy,1:sx)/q(1));
    for k=1:nq
      tq(1:sy,(sx+1):2*sx) = round(t(1:sy,(sx+1):2*sx)/q(3*k-1));
      tq((sy+1):2*sy,1:sx) = round(t((sy+1):2*sy,1:sx)/q(3*k));
      tq((sy+1):2*sy,(sx+1):2*sx) = round(t((sy+1):2*sy,(sx+1):2*sx)/q(3*k+1));
      sy=sy*2;
      sx=sx*2;
    end
  end
end

