clear
close all
clc


%Projet Filtrage et Estimateurs

%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct


% ===== Paramètres
s = signal1.fcno01fz;
NFFT = 512;

%Mise en place de la méthode d'addition-recouvrement
s_rebuilt = addition_recouvrement(s, NFFT, @hamming, 0.5);

%Comparaison signal Pur/Reconstruit
signal_comp = s - s_rebuilt;

%affichage
figure;
plot(s);
title("signal pur");

figure;
plot(s_rebuilt);
title('signal reconstruit');

% Affichage de la comparaison entre le signal pur et le signal reconstruit
% figure;
% plot(signal_comp);
% title('Comparaison: Signal Pur - Signal Reconstruit');

% soundsc(s);
% pause(7);
%soundsc(s_rebuilt);
soundsc(signal_comp);