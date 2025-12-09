clear
close all
clc


%% Projet Filtrage et Estimateurs

%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct

%% Paramètres
s = signal1.fcno01fz;
len_s = length(s);
Fe = 8000; %Hz ,fréquence d'echantillonnage
t=(0:len_s-1)/Fe;
RSB = 10; % = 5:5:15;

% Paramètres du filtre
k_0 = 100;
num = zeros(k_0 + 1, 1);
num(1,1) = 1;
num(k_0 + 1, 1) = 1;
den = 1;

% Paramètres pour le spectrogramme
fenetre = hamming(256); % Fenêtre de Hamming
noverlap = 128; %Chevauchement
nfft = 512; %Nombre de points FFT


%Génération d'un bruit blanc gaussien centré
Bruit = randn(1,len_s);
len_bruit = length(Bruit);


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

%Signal bruité avec ajustement du RSB
signalBruit_ajuste = s + alpha * Bruit'; 

%Mise en place du filtre
%y = filter(num,den,s);
y = filter(num,den,signalBruit_ajuste);

%% Affichage
% Signal de parole pur
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


%Filtre 
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

%% Ecoute des signaux
soundsc(s);
% pause(7);
% soundsc(signalBruit_ajuste);
pause(7);
soundsc(y);