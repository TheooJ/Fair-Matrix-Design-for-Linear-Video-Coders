function y = huffmanlength(p)

% p is given probability vector
%
% huffmanlength(p) is the expected codeword length for the Huffman code
%
% Example: huffmanlength([0.25 0.50 0.25]) is 1.5

N=length(p);
q=sort(p);
S=0;

while N>2
    S=S+q(1)+q(2);
    q=[q(1)+q(2) q(3:N)];
    q=sort(q);
    N=N-1;
end

y=S+1;
