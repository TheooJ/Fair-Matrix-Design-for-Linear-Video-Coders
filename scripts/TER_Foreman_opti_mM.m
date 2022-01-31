%% Vieux code fait pendant le TER dans le cas où on voulait comparer l'initialisation avec les matrices mM
%Utilisé dans comparaison
%% TER - Designing precoding and decoding matrixes
% G : Precoding matrix
% H_i : Decoding matrix at the ith receptor
% Lambda : Source covariance matrix E[tt']
% N : Noise covariance matrix E[vv']

function [G_n, H_n, G_0, H_0, N, MSE] = TER_Mat_Foreman_opti(Lambda,Pe)

load N1
load N2

%% Definition of variables
Lambda = diag(Lambda);
n_ck = length(Lambda); %taille du message a envoyer (et du message estime)
n_se = n_ck; %taille du message une fois precode
n_sr = n_ck; %taille du message recu 

%Pe=2; %Transmitter power budget


%Bruit du recepteur 1
idx_sort=sort(N_canal1)
idx = find(idx_sort<inf);
N_canal1_tronq = idx_sort(idx(1:1+n_ck-1)); %1000
N_1 = diag(N_canal1_tronq);

      %N_1=diag([1 1 2 2 3 3 4 4 5 5 1 1 2 2 3 3 4 4 5 5 1 1 2 2])*10^-2; %Marche bien *1e-7 Sans rien pour droites qui se croisent, 10^-1 pour courbe lisse, 10^-2 et sedumi pour 2 paliers                                                                   
                                                                    
%N_1=diag([1 1 2 2 3 3 4 4 5 5 1 1 2 2 3 3 4 4 5 5 1 2 3 4])*10^-4;

%Bruit du recepteur 2
idx_sort=sort(N_canal1) %1 ou 2
idx = find(idx_sort<inf);
N_canal2_tronq = idx_sort(idx(1:1+n_ck-1)); %3150
N_2 = diag(N_canal2_tronq);

%N_2=diag([1 2 3 4 5 1 2 3 4 5 4 3 2 1 2 3 4 5 4 3 2 1 2 3])*1e-5;% marche bien

%N_2=diag([1 2 3 4 5 1 2 3 4 5 4 3 2 1 2 3 4 5 4 3 2 1 2 3])*10^-4;

%N_2=diag([1 1 2 2 3 3 4 4 5 5 1 1 2 2 3 3 4 4 5 5 1 2 3 4])*10^-3; %Meme que la 1

%N_1 = rand(n_ck)*1e-2; %Ne fonctionne pas, forcement diagonal
%N_2 = rand(n_ck)*1e-2; %Problem with CVX (not solved)

%N_1=diag([1 1 2 2 3 3 4 4 5 5 1 1 2 2 3 3 4 4 5 5 1 2 3 4])*1e-5; %se
%croise

      %N_2=diag([1 2 3 4 5 1 2 3 4 5 4 3 2 1 2 3 4 5 4 3 2 1 2 3])*10^-2; %*1e-7

%N_1=eye(n_ck)*1e-3; %Forcement diagonal *1e-4
%N_2=eye(n_ck)*1e-3; %*1e-4

% N_1= diag(rand(n_ck,1))*1e-3;
% N_2= diag(rand(n_ck,1))*1e-5;

N={N_1  N_2};
MSE = cell(2,1);

%% Initialization
G_n= sqrt(Pe)*eye(n_se,n_ck)/sqrt(trace(Lambda)+1); %+1 %G respecting the power contraint
H_n = cell(1,2);
%W_1 = diag(rand(n_sr,1))*1e-3 + eye(n_sr) ; %rand met en couleur
%W_2 = diag(rand(n_sr,1))*1e-3 + eye(n_sr) ;
W_1 = eye(n_sr);
W_2 = eye(n_sr); %1.2*
W = {W_1 W_2};

G_0 = G_n;
P(1) = trace(G_0*Lambda*G_0'); %Power with G initialization

%% Compute H0
for i = (1:2)
U{1,i} = W{1,i}*G_0;
M{1,i} = U{1,i}*Lambda;
H_0{1,i}=(M{1,i}*U{1,i}' + N{1,i})\M{1,i};
end

eps = 1; %Difference between the norm of two successive iterations
eps_lim = 1e-6; %Stopping criteria
n = 1; %nb of iterations
n_max = 40; %maximum nb of iterations

%% Beginning of the loop
if (P(1) <= Pe) %if the constraint is respected 
    while ( (n < n_max) && (eps >= eps_lim))
        %% Update H with given G
        for i = (1:2)
             U{1,i} = W{1,i}*G_n;
             M{1,i} = U{1,i}*Lambda;
             H_n{1,i}=(M{1,i}*U{1,i}' + N{1,i})\M{1,i};
        end
        %% Solve the problem (CVX) : update G_n_1 given H_n

        cvx_begin
        %cvx_solver sedumi
        variable t_e 
        variable G_n_1(n_se,n_ck) 
        variable Pi_1(n_ck,n_ck)
        variable Pi_2(n_ck,n_ck)
        variable Omega(n_se,n_se)   

        matrix_schur_1 = ([Pi_1+(H_n{1,1}'*W{1,1}*G_n_1)*Lambda+Lambda*(H_n{1,1}'*W{1,1}*G_n_1)' (H_n{1,1}'*W{1,1}*G_n_1) ;
                          (H_n{1,1}'*W{1,1}*G_n_1)' inv(Lambda)]);
        matrix_schur_2 = ([Pi_2+(H_n{1,2}'*W{1,2}*G_n_1)*Lambda+Lambda*(H_n{1,2}'*W{1,2}*G_n_1)' (H_n{1,2}'*W{1,2}*G_n_1) ; 
                          (H_n{1,2}'*W{1,2}*G_n_1)' inv(Lambda)]);
        matrix_schur_3 = ([Omega G_n_1 ;
                           G_n_1' inv(Lambda)]);

        minimize t_e
        trace(Pi_1)+trace(H_n{1,1}'*N{1,1}*H_n{1,1}) <= t_e;
        trace(Pi_2)+trace(H_n{1,2}'*N{1,2}*H_n{1,2}) <= t_e;
        trace(Omega)<=Pe; %Power Budget

        matrix_schur_1 == semidefinite(2*n_ck); %First constraint (first receptor)
        matrix_schur_2 == semidefinite(2*n_ck); %First constraint (second receptor)
        matrix_schur_3 == semidefinite(n_ck+n_se); %Second constraint

        cvx_end;
        
        if (strcmp(cvx_status, ['Solved']) ~= 1) %if status is not 'Solved'
            disp('!!! Problem with CVX');
            pause;
        end
        
        %% Test and update of G_n
        n = n+1;
        epsilon(n) = norm(G_n_1 - G_n,1); %Norm 1
        eps = epsilon(n);
        G_n = G_n_1;
        P(n) = trace(G_n*Lambda*G_n');
        MSE{1}(n) = trace(H_n{1,1}'*W{1,1}*G_n*Lambda*G_n'*W{1,1}'*H_n{1,1}-...
                 H_n{1,1}'*W{1,1}*G_n*Lambda + H_n{1,1}'*N{1,1}*H_n{1,1} - ...
                 Lambda*G_n'*W{1,1}'*H_n{1,1}+Lambda);
        MSE{2}(n) = trace(H_n{1,2}'*W{1,2}*G_n*Lambda*G_n'*W{1,2}'*H_n{1,2}-...
                 H_n{1,2}'*W{1,2}*G_n*Lambda + H_n{1,2}'*N{1,2}*H_n{1,2} - ...
                 Lambda*G_n'*W{1,2}'*H_n{1,2}+Lambda);
    end
    
else
    disp('!!! Initialization of G does not respect the power constraint');
end

% figure;
% plot((2:n), MSE{1}(2:end), (2:n), MSE{2}(2:end));
% legend('MSE_1', 'MSE_2');
end