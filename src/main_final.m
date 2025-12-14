%% Projet Filtrage et Estimateurs
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

clear
clc

%% Premiminaire 2

%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct

%% Paramètres
s = signal1.fcno01fz;
sigma_carre = 100;
len_s = length(s);
Fe = 8000; %Hz ,fréquence d'echantillonnage
t=(0:len_s-1)/Fe;
RSB = 5; % = 5:5:15;

%Paramètres pour le spectrogramme
fenetre = hamming(256); % Fenêtre de Hamming
noverlap = 128; %Chevauchement
nfft = 512; %Nombre de points FFT


%Génération d'un bruit blanc gaussien centré
Bruit = randn(1,len_s);
len_bruit = length(Bruit);

%% Signal bruité sans ajustement ni maitrise du RSB
signalBruit = s +  sigma_carre *Bruit'; % Ajout du bruit au signal avec transposition


%Détermination de la puissance des signaux
%Puissance signal1
Ps= (1/len_s)*sum(s.^2);
fprintf('Ps = %.3f\n',Ps);

%Puissance du bruit
Pb= (1/len_bruit)*sum(Bruit.^2);
fprintf('Pb = %.3f\n',Pb);

%Ajustement
aj = Ps/(Pb * 10^(RSB/10));
alpha = sqrt(aj);

%% Signal bruité avec ajustement du RSB
signalBruit_ajuste = s + alpha * Bruit'; 

%% Affichage
%Signal de parole pur
figure;
subplot(2,1,1);
plot(t,s);
xlabel('Temps (s)')
ylabel('Amplitude')
title('Représentation du signal de parole : fcno01fz')
grid on

subplot(2,1,2);
spectrogram(s, fenetre, noverlap, nfft, Fe, 'yaxis');
title('Spectrogramme du signal de parole : fcno01fz');
colorbar;

%Signal bruité avec ajustement du RSB
figure;
subplot(2,1,1);
plot(t, signalBruit_ajuste);
xlabel('Temps (s)')
ylabel('Amplitude')
title(['Signal de parole bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
grid on;

%Spectrogramme
subplot(2,1,2);
spectrogram(signalBruit_ajuste, fenetre, noverlap, nfft, Fe, 'yaxis');
title(['Spectrogramme du signal bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
colorbar;

%% Préliminaire 3
clear
clc


%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct

%% Paramètres
s = signal1.fcno01fz;
len_s = length(s);
Fe = 8000; %Hz ,fréquence d'echantillonnage
t=(0:len_s-1)/Fe;
RSB = 10; % = 5:5:15;

% Paramètres du filtre
k_0 = 10;
num = zeros(k_0 + 1, 1);
num(1,1) = 1;
num(k_0 + 1, 1) = 1;
den = 1;

% Paramètres pour le spectrogramme
fenetre = hamming(256); % Fenêtre de Hamming
noverlap = 128; %Chevauchement
nfft = 512; %Nombre de points FFT

%Appel à la fonction qui bruite le signal à un RSB donné
[signalBruit_ajuste, bbgc_ajusted] = bruitageSignal(s, RSB);

%Mise en place du filtre
%y = filter(num,den,s);
y = filter(num,den,signalBruit_ajuste);

%% Affichage
% Signal de parole pur
figure;
subplot(2,1,1);
plot(t,s);
xlabel('Temps (s)')
ylabel('Amplitude')
title('Représentation du signal de parole : fcno01fz')
grid on

subplot(2,1,2);
spectrogram(s, fenetre, noverlap, nfft, Fe, 'yaxis');
title('Spectrogramme du signal de parole : fcno01fz');
colorbar;

%Signal Filtré 
figure;
subplot(2,1,1);
plot(t,y);
xlabel('Temps (s)')
ylabel('Amplitude')
title(['Représentation du signal de parole : fcno01fz filtré (k_0 = ' num2str(k_0) ')']);
grid on

subplot(2,1,2);
spectrogram(y,fenetre,noverlap,nfft, Fe, 'yaxis');
title(['Spectrogramme du signal de parole : fcno01fz filtré (k_0 = ' num2str(k_0) ')'])
colorbar;

%% Préliminaire 4
clear
clc

%% Projet Filtrage et Estimateurs

%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct


%% Paramètres
s = signal1.fcno01fz;
len_s = length(s);
Fe = 8000; %Hz ,fréquence d'echantillonnage
t=(0:len_s-1)/Fe;
NFFT = 512;
% type_fenetre : @hamming, @hann, @blackman
type_fenetre = @hann;
overlap = 0.5;


%% Mise en place de la méthode d'addition-recouvrement
s_rebuilt = addition_recouvrement(s, NFFT, type_fenetre, overlap);

%Comparaison signal Pur/Reconstruitx
signal_comp = s - s_rebuilt;

%% Affichage

%Signal original
figure;
subplot(311);
plot(t, s);
xlabel('Temps (s)');
ylabel('Amplitude');
title("Représentation du signal de parole : fcno01fz");
grid on;

%Signal reconstruit
subplot(312);
plot(t, s_rebuilt);
xlabel('Temps (s)');
ylabel('Amplitude');
title(['Représentation du signal de parole reconstruit (Fentrage = Hann et Overlap = ' num2str(overlap) ' )']);
grid on;

%Affichage de la comparaison entre le signal pur et le signal reconstruit
subplot(313);
plot(t, signal_comp);
ylim([-30000, 10000]);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Comparaison entre le signal pur et le signal reconstruit');
grid on;

%% Cahier des charges

%% 2.6 - Base de données de signaux bruités
clear
clc

%Chargement des signaux de parole
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

%% Paramètres de la Soustraction Spectrale
NFFT = 512;
overlap = 0.5;
type_fenetre = @hann;
Fe = 8000;
indice_signal_etude = 3; 


signal_original_1 = x1;
signal_bruite_1 = data_noise(:,indice_signal_etude); % Signal bruité (RSB = 5)
bruit_1 = associated_noise(:, indice_signal_etude); % Bruit associé (RSB = 5)


frames_signal_original_1 = decoupage_trames(signal_original_1, NFFT, overlap);

%Soustraction_spectrale
[signal_reconstruit, frames_signal_bruite , frames_trames_rehausses] = soustraction_spectrale_reconstruction(signal_bruite_1, bruit_1, NFFT, overlap, type_fenetre);

%% Affichages

%Paramètre representation temporelle trame
indice_trame = 1;

%Affichage representation temporelle trame
figure;
subplot(2,1,1);
hold on;
plot(frames_signal_original_1(:,indice_trame), 'b', 'LineWidth', 1.5); %Bleu pour Original
plot(frames_signal_bruite(:,indice_trame), 'r--', 'LineWidth', 1); %Rouge pointillé pour Bruité
plot(frames_trames_rehausses(:,indice_trame), 'g', 'LineWidth', 1.5); %Vert pour Rehaussé
hold off;

title(['Comparaison Temporelle de la Trame ', num2str(indice_trame)]);
xlabel('Échantillons');
ylabel('Amplitude');
legend('Original', 'Bruité', 'Rehaussé');
grid on;

%Paramètre spectres d'amplitudes trames
f = (0:NFFT-1) * (Fe / NFFT);
window = type_fenetre(NFFT);


Amplitude_original = abs(fft(frames_signal_original_1(:, indice_trame) .* window));
Amplitude_bruite = abs(fft(frames_signal_bruite(:, indice_trame) .* window));
Amplitude_rehaussee = abs(fft(frames_trames_rehausses(:, indice_trame) .* window));

%Affichage spectres d'amplitudes trames
subplot(2,1,2);
hold on;
plot(f, Amplitude_original, 'b', 'LineWidth', 1.5);
plot(f, Amplitude_bruite, 'r--', 'LineWidth', 1);
plot(f, Amplitude_rehaussee, 'g', 'LineWidth', 1.5);
hold off;

title(['Comparaison du Spectre d''Amplitude de la Trame ', num2str(indice_trame)]);
xlabel('Fréquence (Hz)');
ylabel('|FFT|');
legend('Original', 'Bruité', 'Rehaussé');
grid on;
axis tight;

%Paramètres affichage signal reconstruit
L=length(signal_reconstruit);
t=(0:L-1)/Fe;
noverlap = 128;

%Signal reconstruit
figure;
subplot(2,1,1);
plot(t, signal_reconstruit);
title('Signal de parole Rehaussé (Soustraction Spectrale)');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on

%Spectrogramme signal reconstruit
subplot(2,1,2);
spectrogram(signal_reconstruit, hamming(256), noverlap, NFFT, Fe, 'yaxis');
title('Spectrogramme du signal de parole rehaussé');
colorbar;

%% 2.6 -Gain en RSB sur fcno01fz
clear
clc


%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat');
x1 = signal1.fcno01fz;  
L = length(x1);

%Paramètres de Soustraction Spectrale
NFFT = 512;
overlap = 0.5;
type_fenetre = @hann;
Fe = 8000;

%Paramètres
N_realisations = 30; % Nombre de réalisations pour la moyenne
RSB_cibles = [5, 10, 15]; 
Nb_RSB = length(RSB_cibles);
D = floor(NFFT/2);

Somme_Gain_RSB = zeros(Nb_RSB, 1); %Matrice pour les resultats

% Détermination des RSB
for r = 1:N_realisations
    for j = 1:Nb_RSB
        
        RSB_cible = RSB_cibles(j);
        
        %Bruitage du signal
        [signal_bruite, bruit_associe] = bruitageSignal(x1, RSB_cible);
        
        %Soustraction Spectrale
        [signal_reconstruit, ~, ~] = soustraction_spectrale_reconstruction(signal_bruite, bruit_associe, NFFT, overlap, type_fenetre); %tilde parce que des paramètres ne nous interessent pas
        
        %Alignement nécessaire
        L_alignee = L-D; 
        s1_aligne = x1((D+1):L); 
        signal_reconstruit_aligne = signal_reconstruit(1:L_alignee); 
        bruit_associe_aligne = bruit_associe((D+1):L);

        %Determination RDB initial
        P_s_initial = mean(s1_aligne.^2);
        P_b_initial = mean(bruit_associe_aligne.^2);
        RSB_initial_dB = 10 * log10(P_s_initial/P_b_initial);
        %Determination RDB final
        bruit_residuel = signal_reconstruit_aligne - s1_aligne; 
        P_b_residuel = mean(bruit_residuel.^2);
        RSB_final_dB = 10 * log10(P_s_initial/P_b_residuel);
        
        %Gain
        Gain_RSB = RSB_final_dB - RSB_initial_dB;
        Somme_Gain_RSB(j) = Somme_Gain_RSB(j)+ Gain_RSB;

    end
end

%Moyenne Finale
Gain_RSB_Moyen = Somme_Gain_RSB/N_realisations;

%% Affichage 
fprintf('Tableau Récapitulatif des Gains en RSB\n');
fprintf(' RSB Initial (dB) | Gain Moyen en RSB (dB)\n');
fprintf('-----------------------------------------\n');
for j = 1:Nb_RSB
    fprintf('%18d | %20.2f\n', RSB_cibles(j), Gain_RSB_Moyen(j));
end
fprintf('-----------------------------------------\n');
