%% Programme global pour rendre la lecture plus facile et
% garder les mêmes composantes de bruit

%% Initialisation

clear all
close all
addpath('wavelets')
addpath('YUV')

N_1=0;
N_2=0;

N={N_1  N_2};

%% Choix de P (seule chose ‡ modifier)
      
%P = [1, 0.1, 0.01, 0.001, 0.0001];
P = [1, 0.1, 0.001];

%Initialisation des PSNR channel 1
A_afficher_ch1_mM = zeros(length(P),1); %sol mM
A_afficher_ch1_ref_1 = zeros(length(P),1); %sol reference 1
A_afficher_ch1_ref_2 = zeros(length(P),1); %sol reference 2

%Initialisation des PSNR channel 2
A_afficher_ch2_mM = zeros(length(P),1);
A_afficher_ch2_ref_1 = zeros(length(P),1);
A_afficher_ch2_ref_2 = zeros(length(P),1);


for j=1:length(P)
    %Calcul des PSNR
    [PSNR_p_1, PSNR_p_2, N] = comparaison_tour(P(j),N,j);
    %Stockage
    A_afficher_ch1_mM(j) = PSNR_p_1{1};
    A_afficher_ch1_ref_1(j) = PSNR_p_1{2};
    A_afficher_ch1_ref_2(j) = PSNR_p_1{3};
    
    A_afficher_ch2_mM(j) = PSNR_p_2{1};
    A_afficher_ch2_ref_1(j) = PSNR_p_2{2};
    A_afficher_ch2_ref_2(j) = PSNR_p_2{3};
end


%Affichage channel 1, min-Max en rouge, lagrangien en bleu

figure(1);
subplot(1,2,1);
semilogx(P(1:length(P)),A_afficher_ch1_mM(1:end),'r*',P(1:length(P)),A_afficher_ch1_ref_1(1:end),'bo',P(1:length(P)),A_afficher_ch1_ref_2(1:end),'bx');
ylabel('PSNR [dB]');
xlabel('Power constraint');
legend('min-Max', 'Reference 1', 'Reference 2');
title('Channel 1');

%Affichage channel 2

subplot(1,2,2);
semilogx(P(1:length(P)),A_afficher_ch2_mM(1:end),'r*',P(1:length(P)),A_afficher_ch2_ref_1(1:end),'bo',P(1:length(P)),A_afficher_ch2_ref_2(1:end),'bx');
ylabel('PSNR [dB]');
xlabel('Power constraint');
legend('min-Max', 'Reference 1', 'Reference 2');
title('Channel 2');
