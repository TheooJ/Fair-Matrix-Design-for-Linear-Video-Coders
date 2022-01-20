%% Fonction d optimisation pour utiliser dans le code global

%% TER - Designing precoding and decoding matrixes with direct optimization

function [G_lagr, H_lagr_1, H_lagr_2, P_1, P_2, Q] = opti_lagr_tour(Lambda,pt,N)

N_1=N{1}; 
N_2=N{2};

%% Definition of variables

n_ck = length(Lambda); %taille du message a envoyer (et du message estime)
n_sc = n_ck;
n_se = n_ck; %taille du message une fois precode
n_sr = n_ck; %taille du message recu

Q=find_perm_inv(Lambda);
Lambda=Q*Lambda*Q';

EQM_1 = cell(1,2);
EQM_2 = cell(1,2);
H_1 = cell(1,2);
H_2 = cell(1,2);
%Perm = cell(1,2);

%G et H renvoyés à la fin ligne 187


%% Find l in reference 1
l1_max=min(n_ck,n_sc) %l_max
l1=0 %Init
m1=[]; %Power allocation for each component of the chunk vector
num1=[]; %Concatenation of following elements for the sum
den1=[];

P_1=find_perm(N_1);
N_1P = P_1*N_1*P_1';
N_2P = P_1*N_2*P_1';

for l=l1_max:-1:1
    sqrt_gamma1 = sum(sqrt(diag(Lambda(1:l,1:l)).*diag(N_1P(1:l,1:l))))...
        / (pt + sum(diag(N_1P(1:l,1:l))));
    m1 = sqrt(diag(Lambda(1:l,1:l)).*diag(N_1P(1:l,1:l)))/sqrt_gamma1 - diag(N_1P(1:l,1:l));
    
    if all(m1>0)
        l1 = length(m1);
        break;
    end
 end

% for i=1:l1_max
%   num1=[num1, sqrt(Lambda(i,i)*N_1(i,i))];
%   sum_num1=sum(num1);
%   den1=[den1, N_1(i,i)];
%   sum_den1=sum(den1);
%   gamma1=((sum_num1)/(pt + sum_den1))^2;
%   m1(i)= sqrt(Lambda(i,i)*N_1(i,i)/gamma1)-N_1(i,i);
%   if (m1(i)>0);
%       l1=l1+1;
%   else
%       break
%   end
% end

%% Compute g1_i for channel and source 1 (reference 1)
g1=[]
for i=1:l1
    g1(i)=sqrt(m1(i))/sqrt(Lambda(i,i));
end

%% Compute the optimal precoding matrix for channel and source 1 
%G_1P=[diag(g1), zeros(l1,n_ck-l1); zeros(n_sc-l1,l1), zeros(n_sc-l1,n_ck-l1)];
%G_1 = G_1P*P_1';

%Dans les autre codes, il y avait G_1=G_1P*P_1, or j'ai déjà pris en compte
%que G_effective=G*P' dans le code du softcast

G_1 = [diag(g1), zeros(l1,n_ck-l1); zeros(n_sc-l1,l1), zeros(n_sc-l1,n_ck-l1)];

%% Compute the optimal decoding matrix for 1, and the decoding matrix for 2
%Only one G and two Hs
H_1{1}=Lambda*G_1'/(G_1*Lambda*G_1' + N_1P); 
H_1{2}=Lambda*G_1'/(G_1*Lambda*G_1' + N_2P);

%% Compute the MSE for ref 1

envoye_ch1=0;
non_envoye_ch1=0;

for i=1:l1
    envoye_ch1=envoye_ch1+sqrt(Lambda(i,i)*N_1P(i,i))
end

for i=l1+1:n_ck
    non_envoye_ch1=non_envoye_ch1+Lambda(i,i)
end

envoye_ch2=0;
non_envoye_ch2=0;

for i=1:l1
    envoye_ch2=envoye_ch2+sqrt(Lambda(i,i)*N_2P(i,i))
end

for i=l1+1:n_ck
    non_envoye_ch2=non_envoye_ch2+Lambda(i,i)
end


EQM_1{1,1} = non_envoye_ch1 + sqrt(sqrt_gamma1^2)*envoye_ch1;  %For channel one
EQM_1{1,2} = non_envoye_ch2 + sqrt(sqrt_gamma1^2)*envoye_ch2;  %For channel two


%% Find l in reference 2
l2_max=min(n_ck,n_sc) %l_max
l2=0 %Init
m2=[]; %Power allocation for each component of the chunk vector
num2=[]; %Concatenation of following elements for the sum
den2=[];

P_2=find_perm(N_2);
N_1P = P_2*N_1*P_2';
N_2P = P_2*N_2*P_2';

for l=l2_max:-1:1
    sqrt_gamma2 = sum(sqrt(diag(Lambda(1:l,1:l)).*diag(N_2P(1:l,1:l))))...
        / (pt + sum(diag(N_2P(1:l,1:l))));
    m2 = sqrt(diag(Lambda(1:l,1:l)).*diag(N_2P(1:l,1:l)))/sqrt_gamma2 - diag(N_2P(1:l,1:l));
    
    if all(m2>0)
        l2 = length(m2);
        break;
    end
end

%% Compute g2_i for channel and source 2 (reference 2)
g2=[]
for i=1:l2
    g2(i)=sqrt(m2(i))/sqrt(Lambda(i,i));
end

%% Compute the optimal precoding matrix for channel and source 2 
%G_2P=[diag(g2), zeros(l2,n_ck-l2); zeros(n_sc-l2,l2), zeros(n_sc-l2,n_ck-l2)];
%G_2 = G_2P*P_2';

G_2 = [diag(g2), zeros(l2,n_ck-l2); zeros(n_sc-l2,l2), zeros(n_sc-l2,n_ck-l2)];

%% Compute the optimal decoding matrix for 2, and the decoding matrix for 1
%Only one G and two Hs
H_2{1}=Lambda*G_2'/(G_2*Lambda*G_2' + N_1P);
H_2{2}=Lambda*G_2'/(G_2*Lambda*G_2' + N_2P);

%% Compute the MSE for ref 2

envoye_ch1=0;
non_envoye_ch1=0;

for i=1:l2
    envoye_ch1=envoye_ch1+sqrt(Lambda(i,i)*N_1P(i,i))
end

for i=l2+1:n_ck
    non_envoye_ch1=non_envoye_ch1+Lambda(i,i)
end

envoye_ch2=0;
non_envoye_ch2=0;

for i=1:l2
    envoye_ch2=envoye_ch2+sqrt(Lambda(i,i)*N_2P(i,i))
end

for i=l2+1:n_ck
    non_envoye_ch2=non_envoye_ch2+Lambda(i,i)
end


EQM_2{1,1} = non_envoye_ch1 + sqrt(sqrt_gamma2^2)*envoye_ch1;  %For channel one
EQM_2{1,2} = non_envoye_ch2 + sqrt(sqrt_gamma2^2)*envoye_ch2;  %For channel two


G_lagr = {G_1 G_2}; 
H_lagr_1 = {H_1{1} H_1{2}}; 
H_lagr_2 = {H_2{1} H_2{2}};
end
