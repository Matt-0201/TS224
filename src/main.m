%% Préparation répertoire de travail
clear 
close all
clc

%% Préliminaire 1
% Paramètres
N = 2048;
sig2 = 1;

% Genration d'un bbgc
bruit = randn(1, N) * sqrt(sig2);

% Fonction d'autocorrélation et estimateurs
xcorr_normalized = xcorr(bruit, N, 'normalized');
xcorr_biased = xcorr(bruit, N, 'biased');
xcorr_unbiased = xcorr(bruit, N, 'unbiased');
freq_corr = -1:1/N:1;

% DSP et spectre de puissance
spectre_puiss = (abs(fftshift(fft(bruit,N))).^2)/N;
DSP = sig2*ones(1,N);
freq_puiss = -1/2:1/N:1/2-1/N;

% Tracés des fonctions d'autocorrélations
figure;
hold on;
plot(freq_corr,xcorr_normalized, 'r', LineWidth=4);
plot(freq_corr,xcorr_biased, 'g', LineWidth=2);
plot(freq_corr,xcorr_unbiased, 'b');
hold off;
legend('Théorique','Biaisé', 'Non biaisé');
title("Tracés de la fonctions d'autocorrélation et des etimateurs d'un BBGC");
ylabel('Amplitutude');
xlabel('Fréquence normalisée')

% Tracés DSP et spectre de puissance
figure;
hold on;
plot(freq_puiss, spectre_puiss, 'g');
plot(freq_puiss, DSP, 'r', LineWidth=2);
hold off;
legend('Spectre de puissance', 'DSP');
title('Tracés de la DSP et du spectre de puissance du BBGC');
ylabel('Amplitutude');
xlabel('Fréquence normalisée')

% Fenetrer: garder les valeurs fiables 

% Périodogramme de Barlett
NFFT = 256;
Fe = 8000;
f = -Fe/2:Fe/NFFT:(Fe/2-Fe/NFFT);

DSP_Barlett = bartlett(bruit, N);

%figure;
%plot(freq_puiss, DSP_Barlett);

% Périodogramme de Welch
DSP_Welch = welch(bruit, NFFT, Fe);

figure;
hold on;
plot(f, DSP_Welch, 'r');
%plot(f, spectre_puiss, 'b');
hold off;

% Corrélogramme
%correlogram = abs(fftshift(fft(xcorr_biased)));

figure;
hold on;
%plot(freq_puiss, correlogram, 'r');
%plot(freq_puiss, spectre_puiss, 'b');
hold off;

% Platitude spectrale
plat_spectrale = geomean(spectre_puiss)/mean(spectre_puiss);


