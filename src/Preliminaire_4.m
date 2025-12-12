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
NFFT = 512;
% type_fenetre : @hamming, @hann, @blackman
type_fenetre = @hann;
overlap = 0.5;


%% Mise en place de la méthode d'addition-recouvrement
s_rebuilt = addition_recouvrement(s, NFFT, type_fenetre, overlap);

%Comparaison signal Pur/Reconstruitx
signal_comp = s - s_rebuilt;

%% Affichage
figure;
subplot(311);
plot(t, s);
xlabel('Temps (s)');
ylabel('Amplitude');
title("Représentation du signal de parole : fcno01fz");
grid on;

subplot(312);
plot(t, s_rebuilt);
xlabel('Temps (s)');
ylabel('Amplitude');
title(['Représentation du signal de parole reconstruit (Fentrage = Hann et Overlap = ' num2str(overlap) ' )']);
grid on;

% Affichage de la comparaison entre le signal pur et le signal reconstruit
subplot(313);
plot(t, signal_comp);
ylim([-30000, 10000]);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Comparaison entre le signal pur et le signal reconstruit');
grid on;

%% Ecoute des signaux de parole
% soundsc(s);
% pause(7);
% soundsc(s_rebuilt);
% pause(7);
% soundsc(signal_comp);