%% Préparation répertoire de travail
clear 
close all
clc

%load("../data/signal_radar_config1.mat")
%addpath("../data");

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

overlap = 0.5;
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

M = 20;                          % Taille de la fenetre de lissage
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

% Platitude spectrale: caractériser le signal

% Expliquer le choix de l'addition recouvrement

%% 2.6 - Base de données de signaux bruités
clc
clear
load("../data/fcno01fz.mat");
load("../data/fcno03fz.mat");
load("../data/fcno04fz.mat");

% Création de la base de donnée de signaux bruités
x1 = fcno01fz;
x2 = fcno03fz;
x3 = fcno04fz;
N_ = length(x3);
x1 = [x1; zeros(N_-length(x1), 1)];  % 0-padding pour mettre les signaux dans une matrice
x2 = [x2; zeros(N_-length(x2), 1)];  % et les filtrer
 
data = [x1, x2, x3];
data_noise = zeros(N_, 9);
associated_noise = zeros(N_, 9);
RSB = 5;
k = 1;
for i=1:9
    [data_noise(:, i), associated_noise(:,i)] = bruitageSignal(data(:,k), RSB);
    RSB = RSB + 5;
    if (RSB == 20)
        RSB = 5;    % On repasse le RSB à 5
        k = k + 1;  % On passe au signal suivant
    end
end

% Méthode de soustraction spectrale
fs = 8000;
DSP_data_noise = abs(fftshift(fft(data_noise(:,1))).^2)/N_;
DSP_noise = abs(fftshift(fft(associated_noise(:,1))).^2)/N_;

DSP_th = abs(fftshift(fft(x1)).^2)/N_;

%stft(data_noise(:,1), fs);

DSP_reconstruct_signal = DSP_data_noise - DSP_noise;

% for i=1:N_
%     new_val = DSP_data_noise(i) - DSP_noise(i);
%     if (new_val < 0)
%         DSP_reconstruct_signal(i) = 0;
%     else
%         DSP_reconstruct_signal(i) = new_val;
%     end
% end

figure;
hold on;
plot(DSP_data_noise,"r", DisplayName="Spectre de puissance signal bruité", LineWidth=2);
plot(DSP_th, "g", DisplayName="SPectre de puissance théorique", LineWidth=3);
plot(DSP_reconstruct_signal, "b", DisplayName="Spectre de puissance du signal reconstruit");
legend();
hold off;

% Affichage des signaux à bruiter (pour vérification)
% figure;
% sgtitle("Premier signal de parole avec bruitage de RSB 5, 10 et 15", fontsize=16);
% subplot(221);
% plot(x1);
% title("Signal de base");
% subplot(222);
% plot(data_noise(:,1));
% title("RSB = 5");
% subplot(223);
% plot(data_noise(:,2));
% title("RSB = 10");
% subplot(224);
% plot(data_noise(:,3));
% title("RSB = 15");
% 
% figure;
% sgtitle("Deuxième signal de parole avec bruitage de RSB 5, 10 et 15", fontsize=16);
% subplot(221);
% plot(x2);
% title("Signal de base");
% subplot(222);
% plot(data_noise(:,4));
% title("RSB = 5");
% subplot(223);
% plot(data_noise(:,5));
% title("RSB = 10");
% subplot(224);
% plot(data_noise(:,6));
% title("RSB = 15");
% 
% figure;
% sgtitle("Troisième signal de parole avec bruitage de RSB 5, 10 et 15", fontsize=16);
% subplot(221);
% plot(x3);
% title("Signal de base");
% subplot(222);
% plot(data_noise(:,7));
% title("RSB = 5");
% subplot(223);
% plot(data_noise(:,8));
% title("RSB = 10");
% subplot(224);
% plot(data_noise(:,8));
% title("RSB = 15");



