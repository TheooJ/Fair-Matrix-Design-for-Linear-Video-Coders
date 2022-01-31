%% Fonction d optimisation pour utiliser dans le code global

%% TER - Designing precoding and decoding matrixes with min-Max solution

function [G_n, H_n] = opti_mM_tour(Lambda,Pe,N,G_init)
        %function [G_n, H_n] = opti_mM_tour(Lambda,Pe,N,G_init)

% N_1=N{1};
% N_2=N{2};

%% Definition of variables

n_ck = length(Lambda); %taille du message a envoyer (et du message estime)
n_se = n_ck; %taille du message une fois precode
n_sr = n_ck; %taille du message recu 

MSE = cell(2,1);

%% Initialization
    G_n= sqrt(Pe)*eye(n_se,n_ck)/sqrt(trace(Lambda)+1); %+1 %G respecting the power contraint
    %G_n=G_init{1}
    
H_n = cell(1,2);

W_1 = eye(n_sr);
W_2 = eye(n_sr); 
W = {W_1 W_2};

G_0 = G_n;
P(1) = trace(G_0*Lambda*G_0'); %Power with G initialization

%% Compute H0
for i = (1:2)
U{1,i} = W{1,i}*G_0;
M{1,i} = U{1,i}*Lambda;
H_0{1,i}=(M{1,i}*U{1,i}' + N{1,i})\M{1,i}; %H!!! et on veut H'. Ici, A\b = inv(A)*b 
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
