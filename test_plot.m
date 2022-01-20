% A_afficher_ch1_mM = A_afficher_ch1_mM';
% A_afficher_ch1_ref_1 = A_afficher_ch1_ref_1';
% A_afficher_ch1_ref_2 = A_afficher_ch1_ref_2';
% 
% A_afficher_ch2_mM = A_afficher_ch2_mM';
% A_afficher_ch2_ref_1 = A_afficher_ch2_ref_1';
% A_afficher_ch2_ref_2 = A_afficher_ch2_ref_2';

%semilogx(P(1:length(P)),A_afficher_ch1_mM(1:end),'r*',P(1:length(P)),A_afficher_ch1_ref_1(1:end),'bo');


%Affichage channel 1

% figure(1);
% subplot(1,2,1);
% semilogx(P(1:length(P)),A_afficher_ch1_mM(1:end),'r*',P(1:length(P)),A_afficher_ch1_ref_1(1:end),'bo',P(1:length(P)),A_afficher_ch1_ref_2(1:end),'bx');
% ylabel('PSNR [dB]');
% xlabel('Power constraint');
% title('Channel 1, mM en rouge, Lagr en bleu');
% 
% %Affichage channel 2
% 
% subplot(1,2,2);
% semilogx(P(1:length(P)),A_afficher_ch2_mM(1:end),'r*',P(1:length(P)),A_afficher_ch2_ref_1(1:end),'bo',P(1:length(P)),A_afficher_ch2_ref_2(1:end),'bx');
% ylabel('PSNR [dB]');
% xlabel('Power constraint');
% title('Channel 2, mM en rouge, Lagr en bleu');

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