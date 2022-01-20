%% TER - Designing precoding and decoding matrixes with min-Max optimization, 
%Etude de l'influence de size(N)
% G : Precoding matrix
% H_i : Decoding matrix at the ith receptor
% Lambda : Source covariance matrix E[tt']
% N_i : Noise covariance matrix E[vv']
% W_i : Matrix between transmitter and i-th receiver

close all;
clear all;

%load N1
%load N2
%load Lambda

%% Definition of variables

%Déclaration au cas par cas
 n_ck = 15; %taille du message a envoyer (et du message estime)
 n_se = n_ck; %taille du message une fois precode
 n_sr = n_ck; %taille du message recu

 
               
 
            N_1=diag([1:(10-1)/(n_ck-1):10])*1e-3;
            N_2=diag([3:(8-3)/(n_ck-1):8])*1e-3;
                    
            
            %Version Kieffer           
            %N_1=diag([1 1 8 8 12 12 16 16 20 20])*1e-3;
            %N_2=diag([1 8 12 1 8 16 20 12 16 20])*1e-3;
            
            
            
            %CHANGER P ICI
            Pe=0.01; %Transmitter power budget
            
            
% Lambda = diag((n_ck:-1:1));
            Lambda = (eye(n_ck)).^3;
            
            
            
%N_2=0.7*N_1;


%N_1=diag([1 1 2 2 3 3 4 4 5 5])*1e-2; %Bruits mm ordre de grandeur
%N_2=diag([1 2 3 4 5 1 2 3 4 5])*1e-2; %si ecart augmente, il faut augmenter Pe     
%N_2=diag([1 2 3 4 5 5 4 3 2 1])*1e-2; %marche
 

%  N_1=diag(rand(n_ck,1))*1e-4;
%  N_2=diag(rand(n_ck,1))*1e-4;

% L_1=sort(L1)
% Lambda = diag(L_1(351:400));
% n_ck = length(Lambda); %taille du message a envoyer (et du message estime)
% n_se = n_ck; %taille du message une fois precode
% n_sr = n_ck; %taille du message recu 

                                                                


%Troncature puis estimation de la covariance
% taille_fen = 18;
% K = 666/taille_fen;
% N_1 = zeros(taille_fen,taille_fen);
% for k = (1:K)
%     x_k = N_canal1_tronq((k-1)*taille_fen+1:k*taille_fen);
%     N_1 = N_1 + x_k*x_k';
% end
% N_1 = N_1/K;

%Prise des N premiers éléments non nuls de N_canal1

% bla=sort(N_canal1);
% idx = find(bla<inf);
% N_canal1_tronq = bla(idx(2000:2019));
% N_1 = diag(N_canal1_tronq);


% bla=sort(N_canal2);
% idx = find(N_canal2<inf);
% N_canal2_tronq = bla(idx(2000:2019));
% N_2 = diag(N_canal2_tronq);


% idx = find(N_canal2<inf);
% N_canal2_tronq = N_canal2(idx(1:10));
% N_2 = diag(N_canal2_tronq);


%N_1=eye(10);
%N_2=eye(10);

N={N_1  N_2};

EQM = cell(1,2);
H_n = cell(1,2);
W_1 = eye(n_sr, n_se) ; %+diag(rand(n_sr,1))*1e-2
W_2 = eye(n_sr, n_se) ; %+diag(rand(n_sr,1))*1e-2 
W = {W_1 W_2};
EQM{1}(1)=0;
EQM{2}(1)=0;

%% Initialization
%G_n= rand(n_se,n_ck)/(n_se*n_ck) ; %G respecting the power contraint
%(never used this line)


%Changer G ICI                        
G_n=sqrt(Pe)*eye(n_se,n_ck)/sqrt(trace(Lambda)); %Initialize the matrix
%just to respect the TPC (usually used)
%G_n=0.5*G_2; %Initialize with the matrix obtained from the optimal solution (ref 1 or 2 TPC)

%Utilisation de la décomposition de Cholesky pour G_n

eps = 100000000; %Difference between the norm of two successive iterations
eps_lim = 10^-3; %Stopping criteria
n = 1; %nb of iterations
n_max = 50; %maximum nb of iterations

G_0=G_n;

P=[];
P(1) = trace(G_0*Lambda*G_0');

%% Beginning of the loop
if (P(1) <= Pe) %if the constraint is respected 
    while ( (n < n_max) && (eps >= eps_lim))
        %% Update H with given G
        for i = (1:2)
            U{1,i} = W{1,i}*G_n;
            M{1,i} = U{1,i}*Lambda;
            H_n{1,i} = (M{1,i}*U{1,i}' + N{1,i})\M{1,i};
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
        
        matrix_schur_1 == semidefinite(2*n_ck); %Premiere contrainte
        matrix_schur_2 == semidefinite(2*n_ck); %Premiere contrainte bis (deuxieme canal)
        matrix_schur_3 == semidefinite(n_ck+n_se); %Deuxieme contrainte
        
        cvx_end;
        
        if (strcmp(cvx_status, ['Solved']) ~= 1) %if status is not 'Solved'
            disp('!!! Problem with CVX');
            pause;
        end
        
        %% Test and update of G_n
        %eps = norm(G_n_1 - G_n,1); %Norm 1
        %epsilon(n) = eps;
        G_n = G_n_1;
        P(n+1) = trace(G_n*Lambda*G_n');
        n = n+1;
        EQM{1}(n) = trace(H_n{1,1}'*W{1,1}*G_n*Lambda*G_n'*W{1,1}'*H_n{1,1}-...
            H_n{1,1}'*W{1,1}*G_n*Lambda + H_n{1,1}'*N{1,1}*H_n{1,1} - ...
            Lambda*G_n'*W{1,1}'*H_n{1,1}+Lambda);
        EQM{2}(n) = trace(H_n{1,2}'*W{1,2}*G_n*Lambda*G_n'*W{1,2}'*H_n{1,2}-...
            H_n{1,2}'*W{1,2}*G_n*Lambda + H_n{1,2}'*N{1,2}*H_n{1,2} - ...
            Lambda*G_n'*W{1,2}'*H_n{1,2}+Lambda);
        eps = max(abs(EQM{1}(n)-EQM{1}(n-1)),abs(EQM{2}(n)-EQM{2}(n-1)));
        epsilon(n) = eps;
    end
    disp('!!! Initialization of G does not respect the power constraint');
end

figure(2);
plot((2:length(EQM{1})), EQM{1}(2:end), (2:length(EQM{2})), EQM{2}(2:end));
legend('Channel 1', 'Channel 2');
%title('Comparaison des MSE');
ylabel('MSE');
xlabel('Number of iterations');

figure(3);
plot((1:length(P)), P(1:length(P)));
ylabel('Power at the transmitter');
xlabel('Number of iterations');
