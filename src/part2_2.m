clear
close all
clc

%Projet Filtrage et Estimateurs

%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct

%Ecoute du signal non bruité
%soundsc(signal1.fcno01fz);

% ===== Paramètres
s = signal1.fcno01fz;
sigma_carre = 100;
len_s = length(s);
Fe = 8000; %Hz ,fréqence d'echantillonnage
t=(0:len_s-1)/Fe;
RSB = 10; % = 5:5:15;

% == Paramètres pour le spectrogramme
fenetre = hamming(256); % Fenêtre de Hamming
noverlap = 128; %Chevauchement
nfft = 512; %Nombre de points FFT


%Génération d'un bruit blanc gaussien centré
Bruit = randn(1,len_s);
len_bruit = length(Bruit);

% ===== Signal bruité sans ajustement ni maitrise du RSB
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

% ===== Signal bruité avec ajustement du RSB
signalBruit_ajuste = s + alpha * Bruit'; 

% ===== Affichage
% == Signal de parole pur
figure;
subplot(2,1,1);
plot(t,s);
%xlim([1 len_s]);
xlabel('Temps (s)')
ylabel('Amplitude')
title('Représentation du signal de parole : fcno01fz')
grid on

subplot(2,1,2);
spectrogram(s, fenetre, noverlap, nfft, Fe, 'yaxis');
title('Spectrogramme du signal de parole : fcno01fz');
colorbar;


% == Bruit BBGC
% figure;
% plot(Bruit);
% xlim([1 len_bruit]);
% grid on;

% == Signal bruité sans ajustmement du RSB
% figure;
% subplot(2,1,1);
% plot(t, signalBruit);
% %xlim([1 len_s]);
% xlabel('Temps (s)')
% ylabel('Amplitude')
% title('Signal de parole bruité sans ajustement du RSB');
% grid on;
% 
% subplot(2,1,2);
% spectrogram(signalBruit, fenetre, noverlap, nfft, Fe, 'yaxis');
% title(['Spectrogramme du signal bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
% colorbar;

%soundsc(signalBruit);

% == Signal bruité avec ajustement du RSB
figure;
subplot(2,1,1);
plot(t, signalBruit_ajuste);
%xlim([1 len_s]);
xlabel('Temps (s)')
ylabel('Amplitude')
title(['Signal de parole bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
grid on;

%soundsc(signalBruit_ajuste);

%Spectrogramme
subplot(2,1,2);
spectrogram(signalBruit_ajuste, fenetre, noverlap, nfft, Fe, 'yaxis');
title(['Spectrogramme du signal bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
colorbar;