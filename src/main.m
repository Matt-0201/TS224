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
xlabel('Fréquence normalisée');

% Périodogramme de Barlett
NFFT = 256;
f = -1/2:1/NFFT:1/2-1/NFFT;

DSP_Bartlett = bartlett(bruit, NFFT);

figure;
subplot(211);
hold on;
plot(f, DSP_Bartlett, 'r', DisplayName='Périodogramme de Bartlett (NFFT points)');
plot(f, DSP(1:NFFT), 'g', DisplayName='Densité spectrale de puissance');
hold off;
legend();
title("Périodogramme de Bartlett et densité spectrale de puissance");
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
subplot(212);
plot(freq_puiss, spectre_puiss, 'b', DisplayName="Spectre de puissance");
legend();
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
title("Spectre de puissance")

% Périodogramme de Welch

overlap = 0.1;
DSP_Welch = welch(bruit, NFFT, overlap);

figure;
subplot(211)
hold on;
plot(f, DSP_Welch, 'r', DisplayName="Périodogramme de Welch (NFFT points)");
plot(f, DSP(1:NFFT), 'g', DisplayName="Densité spectrale de puissance");
hold off;
legend();
title("Périodogramme de Welch et densité spectrale de puissance");
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
subplot(212);
plot(freq_puiss, spectre_puiss, 'b', DisplayName="Spectre de puissance (N points)");
legend();
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
title("Spectre de puissance")

% Périodogramme de Daniell

M = 5;                          % Taille de la fenetre de lissage
DSP_Daniell = daniell(bruit, M);

figure;
subplot(211)
hold on;
plot(freq_puiss, DSP_Daniell, 'r', DisplayName="Périodogramme de Daniell (N points)");
plot(freq_puiss, DSP, 'g', DisplayName="Densité spectrale de puissance");
hold off;
legend();
title("Périodogramme de Daniell et densité spectrale de puissance");
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
subplot(212);
plot(freq_puiss, spectre_puiss, 'b', DisplayName="Spectre de puissance (N points)");
legend();
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
title("Spectre de puissance")

% Corrélogramme

M_correl = N/5;                             % Tronquage, réduit la variance
correlogram = correlogram(bruit, M_correl);
L = length(correlogram);
f_correl = -1/2:1/L:1/2-1/L;

figure;
subplot(211)
hold on;
plot(f_correl, correlogram, 'r', DisplayName="Corrélogramme");
plot(f_correl, DSP(1:L), 'g', DisplayName="Densité spectrale de puissance");
hold off;
legend();
title("Corrélogramme et densité spectrale de puissance");
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
subplot(212);
plot(freq_puiss, spectre_puiss, 'b', DisplayName="Spectre de puissance (N points)");
legend();
ylabel('Amplitutude');
xlabel('Fréquence normalisée');
title("Spectre de puissance")

% Platitude spectrale
%plat_spectrale = geomean(spectre_puiss)/mean(spectre_puiss);


