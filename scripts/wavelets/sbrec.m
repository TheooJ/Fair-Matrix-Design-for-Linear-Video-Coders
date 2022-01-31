function t=sbrec(tq, q)
%sbrec - Do inverse quantization of a subband transformed image
%------------------------------------------------------------------------------
%SYNOPSIS	t = sbrec(tq, q)
%		   where
%
%		     tq   : quantized image to reconstruct
%		     q    : Either a scalar, or a vector of quantization
%                           values, one for each subband of t. The length
%                           of q must be 3n+1, for a positive integer n.
%                           Note that you have top keep track of the number
%                           of subbands splits performed on t yourself, since.
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
%
%SEE ALSO	sbquant, dsbt2, idsbt2, sbdivide, sbmerge
%
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2001 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: sbrec.m,v 1.4 2003/11/18 12:52:23 harna Exp $

  
if (nargin ~= 2)
  error('Wrong number of input arguments.')
end

q = q(:);  % Just for convenience

if iscell(tq)
  t=tq; % Good or bad?
  if length(q)==1
    q=q*ones(length(tq));
  elseif length(q)~=length(tq)
    error('Wrong number of quantization values.')
  end
  for k=1:length(t)
    t{k} = tq{k}*q(k);
  end
else
  nq = floor(length(q)/3);

  if mod(length(q), 3) ~= 1
    error('Wrong number of quantization values.')
  end

  if length(q) == 1
    t=tq*q;   % Only one quantization value
  else
    t=zeros(size(tq));
    sy=floor(size(tq, 1)/(2^nq));
    sx=floor(size(tq, 2)/(2^nq));
    t(1:sy,1:sx) = tq(1:sy,1:sx)*q(1);
    for k=1:nq
      t(1:sy,(sx+1):2*sx) = tq(1:sy,(sx+1):2*sx)*q(3*k-1);
      t((sy+1):2*sy,1:sx) = tq((sy+1):2*sy,1:sx)*q(3*k);
      t((sy+1):2*sy,(sx+1):2*sx) = tq((sy+1):2*sy,(sx+1):2*sx)*q(3*k+1);
      sy=sy*2;
      sx=sx*2;
    end
  end
end

