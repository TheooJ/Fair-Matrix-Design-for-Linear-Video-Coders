function l=arithlength(r)
%arithcode - Arithmetic coding codeword length.
%------------------------------------------------------------------------------
%SYNOPSIS       L = arithlength(r)
%                 r should be a vector describing a histogram of
%                 the symbols of a source. The function estimates the 
%                 number of bits resulting from a memoryless 
%                 arithmetic encoding of the source.
%
%
%SEE ALSO       entropy, huffman, ihist.
%------------------------------------------------------------------------------
%Harald Nautsch                        (C) 2003 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: arithcode.m,v 1.2 2002/12/04 11:51:51 harna Exp harna $

% To simulate the limited precision and rounding effects we'd have in a
% real implementation of an arithmetic coder, we round all probabilities
% to 'prec' bits precision and then assume we need log2(p) bits to code each
% occurence. This is somewhat ad-hoc, but who cares?
  
prec=14;

r=r(find(r));  % Get rid of zero probabilities  
r=sort(r);
ro=r;

r=round(r/sum(r)*2^prec);
r(r==0)=1;     % Make sure we got no new zero probabilities

%Too much trouble. Ignore this, it's just an estimate anyway.
%Adjust the probabilities so they sum to 2^prec
%sumr=sum(r);
%if sumr>2^prec
%  r(end+2^prec-sumr+1:end)=r(end+2^prec-sumr+1:end)-1;
%elseif sumr<2^prec
%  r(end-2^prec+sumr+1:end)=r(end-2^prec+sumr+1:end)+1;
%end

%l=ceil(-sum(ro.*log2(r/2^prec)));

l=ceil(-sum(ro.*log2(r/sum(r))));

