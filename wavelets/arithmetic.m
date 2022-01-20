function C = arithmetic(p,i)
%
% p = given probability vector of length n
%
% For each i between 1 and n, arithmetic(p,i)
% is the i-th arithmetic  codeword
%
n=length(p);
q=cumsum(p);
q=[0 q(1:n-1)];
midpoint=q(i)+p(i)/2;
L=1+ceil(-log2(p(i)));
x=floor(2^L*midpoint);
C=dec2bin(x,L);
