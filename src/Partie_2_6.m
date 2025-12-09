%% 2.6 - Base de données de signaux bruités
clc
clear

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

%% Decoupage en trames 
%On le fait pour le 1er signal, on généralisera apres
s = associated_noise(:,1); %signal d'entrée

%paramètres
len_s = length(s);
Fe = 8000; %Hz ,fréquence d'echantillonnage
t=(0:len_s-1)/Fe;
NFFT = 512;
% type_fenetre : @hamming, @hann, @blackman
type_fenetre = @hann;
overlap = 0.5;


L = length(s);
recouvrement = floor(overlap * NFFT); %Nombre de points qui recouvrent morceaux du signal
ecart = NFFT - recouvrement; %Ecart entre chaque morceau
K = floor((L - NFFT) / ecart) + 1; %nombre total de trames
%s_rebuilt = zeros(L, 1); %signal à reconstruire

fprintf('Nb trames =  %d : NFFT = %d, ecart = %d\n', K, NFFT, ecart);

frames = zeros(NFFT, K);

%% Decoupage en trames 
for i = 1:K %Pour l'instant, on garde la première et la dernière trame
    
    %Découpe en trames
    start_idx = (i-1)*ecart + 1;
    stop_idx = start_idx + NFFT - 1;

    %Extraction + zero padding si nécessaire
    if stop_idx <= L
        frames(:, i) = s(start_idx:stop_idx);
    else %Cas si la dernière trame est incomplete
        trame = zeros(NFFT, 1);
        len_restante = L - start_idx + 1;
        trame(1:len_restante) = s(start_idx:L);
        frames(:, i) = trame;
    end
end
   
% %% Boucle d'Addition - Recouvrement
% for i = 2:K - 1 %On enlève la première et le dernière trame        
%     %Découpe en trames
%     start_idx = (i-1)*ecart + 1;
%     stop_idx = start_idx + NFFT - 1;
% 
%     %Extraction + zero padding si nécessaire
%     if stop_idx <= L
%         trame = s(start_idx:stop_idx);
%     else %Cas si la dernière trame est incomplete
%         trame = zeros(NFFT,1);
%         len_restante = L - start_idx + 1;
%         trame(1:len_restante) = s(start_idx:L); %Le reste de la trame est remplie de 0 (élement neutre)
%     end
% 
%     %Créer une fenêtre de la taille de la trame x
%     window = type_fenetre(length(trame));
% 
%     %Fenetrage de la fenetre
%     trame_windowed = trame .* window;
% 
%     %Reconstitution
%     %s_rebuilt(start_idx:stop_idx) = trame_windowed; 
%     s_rebuilt(start_idx:stop_idx) = s_rebuilt(start_idx:stop_idx) + trame_windowed;
% end

% Calcul de la transformée de Fourier à court terme (STFT)
%[S, F, T] = stft(s, Fe, 'Window', type_fenetre(len_s), 'OverlapLength', overlap*len_s, 'FFTLength', NFFT);



% Méthode de soustraction spectrale
% fs = 8000;
% DSP_data_noise = abs(fftshift(fft(data_noise(:,1))).^2)/N_;
% DSP_noise = abs(fftshift(fft(associated_noise(:,1))).^2)/N_;
% 
% DSP_th = abs(fftshift(fft(x1)).^2)/N_;
% 
% %stft(data_noise(:,1), fs);
% 
% DSP_reconstruct_signal = DSP_data_noise - DSP_noise;

% for i=1:N_
%     new_val = DSP_data_noise(i) - DSP_noise(i);
%     if (new_val < 0)
%         DSP_reconstruct_signal(i) = 0;
%     else
%         DSP_reconstruct_signal(i) = new_val;
%     end
% end

% figure;
% hold on;
% plot(DSP_data_noise,"r", DisplayName="Spectre de puissance signal bruité", LineWidth=2);
% plot(DSP_th, "g", DisplayName="Spectre de puissance théorique", LineWidth=3);
% plot(DSP_reconstruct_signal, "b", DisplayName="Spectre de puissance du signal reconstruit");
% legend();
% hold off;

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



