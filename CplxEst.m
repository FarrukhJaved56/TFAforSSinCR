%% Complexity Analysis
clc
close all
clear all

% Energy Detector NfNC(8NB + 6)
% Auto correlation NfNC(9NB2 +16NB +13/2)
% Spectrogram 4 + Nf (1 + NB +Nfft+NB logNB)+6NBNfft
% Wigner - Ville Nf (2+8Nfft+29NB+NB logNB)-2
% Image processing 14NB2 + 8NB + 51

% Assuming Nf = NB = Nfft
NB = 1:10000;
ED = 8*NB.^2 + 6*NB;
AC = 9*NB.^3 + 16*NB.^2 + 6.5 * NB;
SP = 4 + NB + 8*(NB.^2) + (NB.^2).* log(NB);
WV = 2 * NB + 37 * (NB.^2) + (NB.^2).*log(NB)-2;

axes('fontsize',14,'fontweight','b');
semilogy(NB, AC,'--k',NB, WV,'-.c', NB, SP,':b',NB, ED,'-r','LineWidth',4);
legend('Autocorrelation','WignerVille','Spectrogram','Energy Detector'); 
xlabel('Buffer Size N_B','fontsize',16,'fontweight','b','LineWidth',3); 
ylabel('No of Real Operation','fontsize',16,'fontweight','b','LineWidth',3);
title('COMPLEXITY COMPARISON - SENSING ALGORITHMS','fontsize',18,'fontweight','b','LineWidth',3);
