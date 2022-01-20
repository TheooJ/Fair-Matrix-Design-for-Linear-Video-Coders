%% Comparaison de l'optimisation pour Foreman
% but: comparer mM avec lagrangien pour un cas rŽel

clear all
close all
addpath('wavelets')
addpath('YUV')

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


                    % Change P here
                            P=0.01;
                    

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
    
    %% Min-max vs lagrangian optimization (finding G, Hi and Ni)
    
    
                    
    [G_mM, H_mM, N, MSE] = TER_Foreman_opti_mM_perm(Lambda,P);                                          
    [G_lagr, H_lagr_1, H_lagr_2, MSE_1, MSE_2, P_1, P_2, Q] = TER_Foreman_opti_lagr_new_perm(Lambda,P,N); %G_lagr et G_mM cell array size 2
    
    
    
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
                    
        %s = G*x;
                s = Pio'*G*Qio'*x; %Rajouté

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
                z_1=Pio*y_1; %Rajouté
        
        y_2 = s + b_2;
                z_2=Pio*y_2; %Rajouté

        %% Decoding 
        
        
                         
%         x_hat_1 = H{1}'*y_1; 
%         x_hat_2 = H{2}'*y_2;
                x_hat_1 = Qio*H{1}*z_1; %Rajouté
                x_hat_2 = Qio*H{2}*z_2; 
                 


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
    
        %Channel1
        figure(1);
        subplot(3,1,1);
        imshow(RGB);
        title('Image at the transmitter');
        subplot(3,1,2);
        imshow(RGB_hat_p_1{1});
        title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_1{1})]);
        subplot(3,1,3);
        imshow(RGB_hat_p_1{2});
        title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_1{2})]);
        sgtitle('Reference 1, Channel 1');
        
        
        %Channel2
        figure(2);
        subplot(3,1,1);
        imshow(RGB);
        title('Image at the transmitter');
        subplot(3,1,2);
        imshow(RGB_hat_p_2{1});
        title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_2{1})]);
        subplot(3,1,3);
        imshow(RGB_hat_p_2{2});
        title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_2{2})]);
        sgtitle('Reference 1, Channel 2');
        
        %REFERENCE 2, G et H lagrangien optimisé pour Channel 2
        
        %p_1 ch1 
        %{1} mM 
        %{2} lagrangien ref 1
        %{3} lagrangien ref 2
        
        %Channel1
        figure(3);
        subplot(3,1,1);
        imshow(RGB);
        title('Image at the transmitter');
        subplot(3,1,2);
        imshow(RGB_hat_p_1{1}); %p_1{1} mM ch1
        title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_1{1})]);
        subplot(3,1,3);
        imshow(RGB_hat_p_1{3});
        title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_1{3})]);
        sgtitle('Reference 2, Channel 1');
        
        
        %Channel2
        figure(4);
        subplot(3,1,1);
        imshow(RGB);
        title('Image at the transmitter');
        subplot(3,1,2);
        imshow(RGB_hat_p_2{1}); %p_2{1} mM ch2
        title(['Image estimated with optimized mM matrices, PSNR = ' num2str(PSNR_p_2{1})]);
        subplot(3,1,3);
        imshow(RGB_hat_p_2{3});
        title(['Image estimated with optimized lagrangian matrices, PSNR = ' num2str(PSNR_p_2{3})]);
        sgtitle('Reference 2, Channel 2');
        
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
    
        
        %Plot MSE si même longueur
    
%     figure(3);
%     plot((2:length(MSE{1})), MSE{1}(2:end), (2:length(MSE{2})), MSE{2}(2:end));
%     legend('Channel 1', 'Channel 2');
%     ylabel('MSE');
%     xlabel('Number of iterations');

end

% %We get the PSNR of each image and each channel
% for c = (1:nb_im) %init ch1
%     inter = PSNR_cell_1{c}
%     PSNR_ch_1_init(c) = inter{1}
% end
% 
% for c = (1:nb_im) %init ch2
%     inter = PSNR_cell_2{c}
%     PSNR_ch_2_init(c) = inter{1}
% end
% 
% for c = (1:nb_im) %opti ch1
%     inter = PSNR_cell_1{c}
%     PSNR_ch_1_opti(c) = inter{2}
% end
% 
% for c = (1:nb_im) %opti ch2
%     inter = PSNR_cell_2{c}
%     PSNR_ch_2_opti(c) = inter{2}
% end

% % PSNR ch1
% figure(56);
% plot((1:nb_im), PSNR_ch_1_opti-PSNR_ch_1_init);
% legend('Difference of PSNR (ch1)');
% xlabel('n-th image')
% ylabel('PSNR [dB]');
% 
% % PSNR ch2
% figure(64);
% plot((1:nb_im), PSNR_ch_2_opti-PSNR_ch_2_init);
% legend('Difference of PSNR (ch2)');
% xlabel('n-th image')
% ylabel('PSNR [dB]');
% 
% % Both PSNR
% figure(58);
% plot((1:nb_im), PSNR_ch_1_opti-PSNR_ch_1_init, (1:nb_im), PSNR_ch_2_opti-PSNR_ch_2_init);
% legend('Difference of PSNR (ch1)', 'Difference of PSNR (ch2)');
% xlabel('n-th image');
% ylabel('PSNR [dB]');

% %%  We get the images
% for c = (1:nb_im) %opti ch1
%     inter = RGB_hat_cell_1{c};
%     RGB_ch_1_opti{c} = inter{2};
% end
% 
% for c = (1:nb_im) %opti ch2
%     inter = RGB_hat_cell_2{c};
%     RGB_ch_2_opti{c} = inter{2};
% end
% 
% for c = (1:nb_im) %init ch1
%     inter = RGB_hat_cell_1{c};
%     RGB_ch_1_init{c} = inter{1};
% end
% 
% for c = (1:nb_im) %init ch2
%     inter = RGB_hat_cell_2{c};
%     RGB_ch_2_init{c} = inter{1};
% end

% %% Displaying videos
% for im_idx = 1:nb_im
%     figure(90);
%     subplot(2,2,1);
%     imshow(RGB_ch_1_init{im_idx});
%     title(['Ch1 init, PSNR = ' num2str(PSNR_ch_1_init(im_idx))]);
%     subplot(2,2,2);
%     imshow(RGB_ch_2_init{im_idx});
%     title(['Ch2 init, PSNR = ' num2str(PSNR_ch_2_init(im_idx))]);
%     subplot(2,2,3);
%     imshow(RGB_ch_1_opti{im_idx});
%     title(['Ch1 opti, PSNR = ' num2str(PSNR_ch_1_opti(im_idx))]);
%     subplot(2,2,4);
%     imshow(RGB_ch_2_opti{im_idx});
%     title(['Ch2 opti, PSNR = ' num2str(PSNR_ch_2_opti(im_idx))]);
% end
%% 
fclose(fid);
