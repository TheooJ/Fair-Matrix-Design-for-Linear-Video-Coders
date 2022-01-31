%% Fonction comparaison de l optimisation pour Foreman,
% but: comparer mM avec lagrangien pour un cas rŽel

function [PSNR_p_1, PSNR_p_2 ,N] = comparaison_tour(P,N,K)

load N1;
load N2;

%% Loading data

% Size of the chunks
rck = 36;
cck = 44;
% rck = 24;
% cck = 44;

% Opening the video file
fid = fopen('wavelets/Foreman.qcif');

% Number of images to process
nb_im = 1;

PSNR = cell(3,1);
RGB_hat = cell(3,1);

%% Creation of cells containing data associated with every image
%Msg
% YUV Components of the image
compY_cell = cell(nb_im,1);
compU_cell = cell(nb_im,1);
compV_cell = cell(nb_im,1);

% Mean of the YUV components of the image
m_compY_cell = cell(nb_im,1);
m_compU_cell = cell(nb_im,1);
m_compV_cell = cell(nb_im,1);

% RGB images (at the transmitter)
RGB_cell = cell(nb_im,1);

% DCT2D of the YUV components of the image
dct_compY_cell = cell(nb_im,1);
dct_compU_cell = cell(nb_im,1);
dct_compV_cell = cell(nb_im,1);

% Vector to send and its covariance matrix
x_cell = cell(nb_im,1);
Lambda_cell = cell(nb_im,1);

%Ch1
% Estimated vector
x_hat_cell_1 = cell(nb_im,1);

% Estimated DCT2D of the image
dct_compY_hat_cell_1 = cell(nb_im,1);
dct_compU_hat_cell_1 = cell(nb_im,1);
dct_compV_hat_cell_1 = cell(nb_im,1);

% Estimated YUV components of the image
compY_hat_cell_1 = cell(nb_im,1);
compU_hat_cell_1 = cell(nb_im,1);
compV_hat_cell_1 = cell(nb_im,1);

% PSNR
PSNR_cell_1 = cell(nb_im,1);

% Estimated RGB image
RGB_hat_cell_1 = cell(nb_im,1);

%Ch2
% Estimated vector
x_hat_cell_2 = cell(nb_im,1);

% Estimated DCT2D of the image
dct_compY_hat_cell_2 = cell(nb_im,1);
dct_compU_hat_cell_2 = cell(nb_im,1);
dct_compV_hat_cell_2 = cell(nb_im,1);

% Estimated YUV components of the image
compY_hat_cell_2 = cell(nb_im,1);
compU_hat_cell_2 = cell(nb_im,1);
compV_hat_cell_2 = cell(nb_im,1);

% PSNR
PSNR_cell_2 = cell(nb_im,1);

% Estimated RGB image
RGB_hat_cell_2 = cell(nb_im,1);

%% Debut boucle sur les images
for im_idx = 1:nb_im
    
    %% Mise en forme image
    % Loading the image
    [compY,compU,compV]=yuv_readimage(fid,'qcif');
    
    % Size of the Y component
    [nrows,ncols] = size(compY);
    
    % Making U and V components computable
    compU = compU';
    compV = compV';
    
    % Creating a RGB image
    RGB = yuv2rgb(compY,compU,compV,'YUV420_8');

    % Computing the mean of each component
    m_compY = mean(compY(:));
    m_compU = mean(compU(:));
    m_compV = mean(compV(:));
    
    % Computing the DCT-2D
    dct_compY = dct2(compY-m_compY);
    dct_compU = dct2(compU-m_compU);
    dct_compV = dct2(compV-m_compV);
    
    %% Construction of the vectorized chunks matrix
    x=[];
    for i=1:(ncols/cck)
        for j=1:(nrows/rck)
            X = dct_compY(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck);
            x = [x;X(:)'];
        end
    end
    for i=1:(ncols/2/cck)
        for j=1:(nrows/2/rck)
            X = dct_compU(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck);
            x = [x;X(:)'];
        end
    end
    for i=1:(ncols/2/cck)
        for j=1:(nrows/2/rck)
            X = dct_compV(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck);
            x = [x;X(:)'];
        end
    end

    % Variance of each chunk
    Lambda = var(x');
    
    %% Pour calculer le bruit à la première itération puis garder le même 
    Lambda = diag(Lambda);
    n_ck = length(Lambda); %taille du message a envoyer (et du message estime)
    
    if K==1
        %Bruit du recepteur 1
        N_non_inf = N_canal1(N_canal1<inf);
        
        for i = 1:n_ck
            ind = randperm(numel(N_non_inf), 1); % select one element out of numel(x) elements, with the probability of occurrence of the element in x
            N_vect_1(i) = N_non_inf(ind);
            N_non_inf(ind) = []; % delete this element from the sample, such that the picked elements are unique
        end
        
        N_1 = diag(N_vect_1);
        
        
        %Bruit du recepteur 2
        N_non_inf = N_canal2(N_canal2<inf);
        
        for i = 1:n_ck
            ind = randperm(numel(N_non_inf), 1); % select one element out of numel(x) elements, with the probability of occurrence of the element in x
            N_vect_2(i) = N_non_inf(ind);
            N_non_inf(ind) = []; % delete this element from the sample, such that the picked elements are unique
        end
        
        N_2 = diag(N_vect_2);
                
    else 
        N_1=N{1};
        N_2=N{2};
    end
    
    N={N_1  N_2};

    
    %% Min-max vs lagrangian optimization (finding G, Hi and Ni)
    
    [G_lagr, H_lagr_1, H_lagr_2, P_1, P_2, Q] = opti_lagr_tour(Lambda,P,N); %G_lagr et G_mM cell array size 2
   
    %Initialisation avec G identité pondérée
    %[G_mM, H_mM] = opti_mM_tour(Lambda,P,N);                                         
    
    %Initialisation avec G 
    G_lagr_des={G_lagr{1}*Q G_lagr{2}*Q}
    [G_mM, H_mM] = opti_mM_tour(Lambda,P,N,G_lagr_des); 
    %Choix init ligne 20 dans opti_mM_tour

%             G_1=G_lagr{1};
%             G_2=G_lagr{2};
%             CmM=trace(G_mM*Lambda*G_mM');
%             C1=trace(G_1*Q*Lambda*Q'*G_1');
%             C2=trace(G_2*Q*Lambda*Q'*G_2');
% 
%             if (C1 <= P) && (C2 <= P) && (CmM <=P)
%                 continue 
%             else
%                 disp('Initialization of G does not respect the power constraint')
%             end
%
%             pause;
    
    for l=(1:3) %loop to compare both optimizations
        if l == 1 %Cas minMax
            G = G_mM;
            H{1} = (H_mM{1})'; %On transpose ici et pas dans le lagrangien
            H{2} = (H_mM{2})'; %pour avoir t=H^h*z plus bas 
            Pio=eye(24);
            Qio=eye(24);
        elseif l == 2 %Cas lagrangien optimisé pour channel 1
            G = G_lagr{1};
            H = H_lagr_1;
            Pio=P_1;
            Qio=Q;
        elseif l == 3 %Cas lagrangien optimisé pour channel 2
            G = G_lagr{2};
            H = H_lagr_2;
            Pio=P_2;
            Qio=Q;
        end
        
        %% Precoding
                    
        s = Pio'*G*Qio*x; 

        %% Simulation of the transmission channel
        % Gaussian noise with a given covariance matrix at the first receiver
        b_1 = randn(size(x));
        Q_1 = chol(N{1});
        b_1 = Q_1*b_1;

        % Gaussian noise with a given covariance matrix at the second receiver
        b_2 = randn(size(x));
        Q_2 = chol(N{2});
        b_2 = Q_2*b_2;

        % Received vectors
        y_1 = s + b_1;
        z_1=Pio*y_1; 
        
        y_2 = s + b_2;
        z_2=Pio*y_2; 

        %% Decoding 
        
        x_hat_1 = Qio'*H{1}*z_1;
        x_hat_2 = Qio'*H{2}*z_2;

        %% Reconstruction of the vectorized chunks matrix
        %Channel1
        dct_compY_hat_1 = zeros(nrows,ncols);
        dct_compU_hat_1 = zeros(nrows/2,ncols/2);
        dct_compV_hat_1 = zeros(nrows/2,ncols/2);

        % Luma Y
        for i=1:(ncols/cck)
            for j=1:(nrows/rck)
                idx = j+(i-1)*(nrows/rck);
                X_1 = x_hat_1(idx,:);
                dct_compY_hat_1(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_1,rck,cck);
            end
        end
        % Chroma U
        for i=1:(ncols/2/cck)
            for j=1:(nrows/2/rck)
                idx = j+(i-1)*(nrows/2/rck)+(ncols/cck)*(nrows/rck);
                X_1 = x_hat_1(idx,:);
                dct_compU_hat_1(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_1,rck,cck);
            end
        end
        % Chroma V
        for i=1:(ncols/2/cck)
            for j=1:(nrows/2/rck)
                idx = j+(i-1)*(nrows/2/rck)+(ncols/cck)*(nrows/rck)+(ncols/2/cck)*(nrows/2/rck);
                X_1 = x_hat_1(idx,:);
                dct_compV_hat_1(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_1,rck,cck);
            end
        end

        %Channel2
        dct_compY_hat_2 = zeros(nrows,ncols);
        dct_compU_hat_2 = zeros(nrows/2,ncols/2);
        dct_compV_hat_2 = zeros(nrows/2,ncols/2);

        % Luma Y
        for i=1:(ncols/cck)
            for j=1:(nrows/rck)
                idx = j+(i-1)*(nrows/rck);
                X_2 = x_hat_2(idx,:);
                dct_compY_hat_2(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_2,rck,cck);
            end
        end
        % Chroma U
        for i=1:(ncols/2/cck)
            for j=1:(nrows/2/rck)
                idx = j+(i-1)*(nrows/2/rck)+(ncols/cck)*(nrows/rck);
                X_2 = x_hat_2(idx,:);
                dct_compU_hat_2(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_2,rck,cck);
            end
        end
        % Chroma V
        for i=1:(ncols/2/cck)
            for j=1:(nrows/2/rck)
                idx = j+(i-1)*(nrows/2/rck)+(ncols/cck)*(nrows/rck)+(ncols/2/cck)*(nrows/2/rck);
                X_2 = x_hat_2(idx,:);
                dct_compV_hat_2(1+(j-1)*rck:j*rck,1+(i-1)*cck:i*cck) = reshape(X_2,rck,cck);
            end
        end

        %% Computing of the inverse DCT-2D
        %Channel1
        compY_hat_1 = idct2(dct_compY_hat_1)+m_compY;
        compU_hat_1 = idct2(dct_compU_hat_1)+m_compU;
        compV_hat_1 = idct2(dct_compV_hat_1)+m_compV;

        PSNR_1 = 10*log10(255^2/var(compY_hat_1(:)-compY(:)));
        PSNR_p_1{l} = PSNR_1;

        RGB_hat_1 = yuv2rgb(compY_hat_1,compU_hat_1,compV_hat_1,'YUV420_8');
        RGB_hat_p_1{l} = RGB_hat_1;

        %Channel2
        compY_hat_2 = idct2(dct_compY_hat_2)+m_compY;
        compU_hat_2 = idct2(dct_compU_hat_2)+m_compU;
        compV_hat_2 = idct2(dct_compV_hat_2)+m_compV;

        PSNR_2 = 10*log10(255^2/var(compY_hat_2(:)-compY(:)));
        PSNR_p_2{l} = PSNR_2;

        RGB_hat_2 = yuv2rgb(compY_hat_2,compU_hat_2,compV_hat_2,'YUV420_8');
        RGB_hat_p_2{l} = RGB_hat_2;


    end %end for initi or opti
    
    %REFERENCE 1, G et H lagrangien optimisé pour Channel 1
    
%         %Channel1
%         figure(1);
%         subplot(3,1,1);
%         imshow(RGB);
%         title('Image at the transmitter');
%         subplot(3,1,2);
%         imshow(RGB_hat_p_1{1});
%         title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_1{1})]);
%         subplot(3,1,3);
%         imshow(RGB_hat_p_1{2});
%         title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_1{2})]);
%         sgtitle('Reference 1, Channel 1');
%         
%         
%         %Channel2
%         figure(2);
%         subplot(3,1,1);
%         imshow(RGB);
%         title('Image at the transmitter');
%         subplot(3,1,2);
%         imshow(RGB_hat_p_2{1});
%         title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_2{1})]);
%         subplot(3,1,3);
%         imshow(RGB_hat_p_2{2});
%         title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_2{2})]);
%         sgtitle('Reference 1, Channel 2');
%         
%         %REFERENCE 2, G et H lagrangien optimisé pour Channel 2
%         
%         %p_1 ch1 
%         %{1} mM 
%         %{2} lagrangien ref 1
%         %{3} lagrangien ref 2
%         
%         %Channel1
%         figure(3);
%         subplot(3,1,1);
%         imshow(RGB);
%         title('Image at the transmitter');
%         subplot(3,1,2);
%         imshow(RGB_hat_p_1{1}); %p_1{1} mM ch1
%         title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_1{1})]);
%         subplot(3,1,3);
%         imshow(RGB_hat_p_1{3});
%         title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_1{3})]);
%         sgtitle('Reference 2, Channel 1');
%         
%         
%         %Channel2
%         figure(4);
%         subplot(3,1,1);
%         imshow(RGB);
%         title('Image at the transmitter');
%         subplot(3,1,2);
%         imshow(RGB_hat_p_2{1}); %p_2{1} mM ch2
%         title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_2{1})]);
%         subplot(3,1,3);
%         imshow(RGB_hat_p_2{3});
%         title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_2{3})]);
%         sgtitle('Reference 2, Channel 2');
        
    %% Saving the datas of the image
    %Msg
    compY_cell{im_idx} = compY;
    compU_cell{im_idx} = compU;
    compV_cell{im_idx} = compV;
    
    RGB_cell{im_idx} = RGB;
    
    m_compY_cell{im_idx} = m_compY;
    m_compU_cell{im_idx} = m_compU;
    m_compV_cell{im_idx} = m_compV;
    
    dct_compY_cell{im_idx} = dct_compY;
    dct_compU_cell{im_idx} = dct_compU;
    dct_compV_cell{im_idx} = dct_compV;
    
    x_cell{im_idx} = x;
    Lambda_cell{im_idx} = Lambda;
    
    %Ch1
    x_hat_cell_1{im_idx} = x_hat_1;
    
    dct_compY_hat_cell_1{im_idx} = dct_compY_hat_1;
    dct_compU_hat_cell_1{im_idx} = dct_compU_hat_1;
    dct_compV_hat_cell_1{im_idx} = dct_compV_hat_1;
    
    compY_hat_cell_1{im_idx} = compY_hat_cell_1;
    compU_hat_cell_1{im_idx} = compU_hat_cell_1;
    compV_hat_cell_1{im_idx} = compV_hat_cell_1;
    
    PSNR_cell_1{im_idx} = PSNR_p_1;
        
    RGB_hat_cell_1{im_idx} = RGB_hat_p_1;
    
    %Ch2    
    x_hat_cell_2{im_idx} = x_hat_2;
    
    dct_compY_hat_cell_2{im_idx} = dct_compY_hat_2;
    dct_compU_hat_cell_2{im_idx} = dct_compU_hat_2;
    dct_compV_hat_cell_2{im_idx} = dct_compV_hat_2;
    
    compY_hat_cell_2{im_idx} = compY_hat_cell_2;
    compU_hat_cell_2{im_idx} = compU_hat_cell_2;
    compV_hat_cell_2{im_idx} = compV_hat_cell_2;
    
    PSNR_cell_2{im_idx} = PSNR_p_2;
        
    RGB_hat_cell_2{im_idx} = RGB_hat_p_2;
    

end

end
