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
len_signal1 = length(s);
RSB = 5; % = 5:5:15;

%Génération d'un bruit blanc gaussien centré
Bruit = randn(1,len_signal1);
len_bruit = length(Bruit);

% ===== Signal bruité sans ajustement ni maitrise du RSB
signalBruit = s +  sigma_carre *Bruit'; % Ajout du bruit au signal avec transposition


%Détermination de la puissance des signaux
%Puissance signal1
Ps= (1/len_signal1)*sum(s.^2);
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
%Signal de parole
figure;
plot(s)
title('Représentation du 1er signal de parole')
grid on

%Bruit BBGC
% figure;
% plot(Bruit);
% grid on;

%Signal bruité sans ajustmement du RSB
figure;
plot(signalBruit);
title('Signal de parole bruité sans ajustement du RSB');
grid on;

%soundsc(signalBruit);

%Signal bruité avec ajustement du RSB
figure;
plot(signalBruit_ajuste);
title(['Signal de parole bruité avec ajustement du RSB (RSB = ' num2str(RSB) ' dB)']);
grid on;

%soundsc(signalBruit_ajuste);