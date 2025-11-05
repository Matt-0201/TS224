clear
close all
clc

%Projet Filtrage et Estimateurs
%Chargement d'un signal de parole
signal1 = load('../data/fcno01fz.mat'); %Attention, c'est un struct

%Paramètres
sigma_carre = 0.5;
len_signal1 = length(signal1.fcno01fz)

%Génération d'un bruit blanc gaussien centré
Bruit = sigma_carre * randn(1,len_signal1);

%Signal bruité
signalBruit = signal1.fcno01fz + Bruit; % Ajout du bruit au signal %Il y a un transpose a faire 

%Affichage

%Signal de parole
figure;
plot(signal1.fcno01fz)
title('Représentation du 1er signal de parole')
grid on

%Bruit BBGC
figure;
plot(Bruit);
grid on;

% Signal bruité
figure;
plot(signalBruit);
title('Signal de parole bruité');
grid on;