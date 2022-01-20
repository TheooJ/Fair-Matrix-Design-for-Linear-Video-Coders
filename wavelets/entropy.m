function y = entropy(x)
% The name of this m-file is entropy.m
% x is a data vector
% y=entropy(x) is the first-order entropy of x

P=frequency(x)/length(x);
y=sum(-P.*log2(P));

