function R = runlengths(M)
[s,r]=size(M);
for i=1:s;
x(r*(i-1)+1:r*i)=M(i,:);
end
N=r*s;
y=x(2:N);
x1=x(1:N-1);
z=y+x1;
j=find(z==1);
j1=[j N];
j2=[0 j];
R=j1-j2;



