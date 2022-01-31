%Let  x  be a bitstring and suppose  x  is the i-th string in
%the list of all bitstrings 0, 1, 00, 01, 10, 11, 000, ...
% Then, bitstring_to_index(x) = i.  The integer  i  shall be
%called the index of x.
function y=bitstring_to_index(x)
N=length(x);
S=1;
for i=1:N;
S=2*S+x(i);
end
y=S-1;




