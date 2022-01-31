%% Fonction qui trouve la matrice de permutation inverse

function [P] = find_perm_inv(N); %N matriciel

N_vect=N(:); %On vectorise N
i=(N_vect==0); %On cherche les composants nuls
N_vect(i)=[]; %On les retire
N=N_vect';

P=eye(length(N));

i=1;

while (not(isequal(N,sort(N))) && (i<length(N)))
    if N(i)<N(i+1);
        N([i;i+1])=N([i+1;i]);
        P([i i+1],:)=P([i+1 i],:);
        i=1;
    else i=i+1;
    end
end

end
