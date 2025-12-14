function [signal_rehausse, frames_signal_bruite, trames_restaurees] = soustraction_spectrale_reconstruction(signal_bruite, bruit, NFFT, overlap, type_fenetre)
    % Cette fonction applique la méthode de rehaussement de la parole basé sur la méthode de la soustraction spectrale, suivi d'une reconstruction temporelle par addition–recouvrement.
    
    %Entrées :
    %signal_bruite : signal de parole bruité (vecteur colonne)
    %bruit         : bruit associé au signal de parole
    %NFFT          : taille des trames utilisées pour l'analyse spectrale
    %overlap       : taux de recouvrement entre deux trames (entre 0 et 1)
    %type_fenetre  : type de fenêtre appliquée aux trames (@hamming, @hann, @blackman, etc.)
    
    % Sorties :
    %signal_rehausse        : signal de parole rehaussé après reconstruction
    %frames_signal_bruite   : trames du signal bruité
    %trames_restaurees      : trames temporelles restaurées après soustraction spectrale et transformée de Fourier inverse

    %% Paramètres
    L= length(signal_bruite);
    window = type_fenetre(NFFT);

    %% Decoupage en trames
    frames_signal_bruite = decoupage_trames(signal_bruite, NFFT, overlap);
    frames_bruit = decoupage_trames(bruit, NFFT, overlap);

    K = size(frames_signal_bruite, 2); %Nombre de trames

    %A conserver pour la reconstruction
    fft_frames_signal_bruite = zeros(NFFT, K);
    fft_frames_bruit = zeros(NFFT, K);
    
    %% Transformee de Fourrier de toutes les trames, signal bruité + bruit
    %Comme on a deja des trames, faire la tranformee de fourrier sur chacune des trames, revient a faire un stft
    for i = 1:K
        fft_frames_signal_bruite(:, i) = fft(frames_signal_bruite(:, i) .* window);
        fft_frames_bruit(:, i) = fft(frames_bruit(:, i) .* window);
    end
    
    %% Puissance spectrale
    %Calcul de la puissance spectrale
    P_frames_signal_bruite = abs(fft_frames_signal_bruite).^2;
    P_frames_bruit = abs(fft_frames_bruit).^2;


    %% Soustraction spectrale 
    %Calcul de la soustraction spectrale pour toutes les trames
    P_estim_parole = P_frames_signal_bruite - P_frames_bruit;

    %% Half-Wave Rectification
    %Si P_estim_parole est négatif, on le met à zéro.
    P_estim_parole(P_estim_parole < 0) = 0;

    %Reconstruction du Spectre Complexe Estimé
    %Transformée de Fourrier inverse
    %Amplitude estimée
    Amplitude_estim_parole = sqrt(P_estim_parole);
 
    %Extraction phase du signal bruité
    Phase_bruitee = angle(fft_frames_signal_bruite);
 
    %Reconstruction spectre complexe estimé
    Spectre_complexe_estim_parole = Amplitude_estim_parole .* exp(1i*Phase_bruitee);
 
    %% Transformée de Fourier Inverse sur le spectre complexe corrigé
    %trames_restaurees = ifft(Spectre_complexe_estim_parole, [], 1);
    trames_restaurees = real(ifft(Spectre_complexe_estim_parole, [], 1));

    %% Méthode addition recouvrement
    recouvrement = floor(overlap * NFFT); %Nombre de points qui recouvrent morceau du signal
    ecart = NFFT - recouvrement; %Ecart entre chaque morceau
    K = floor((L - NFFT) / ecart) + 1; %nombre total de trames
    s_rebuilt = zeros(L, 1); %signal à reconstruire
    
    %% Boucle d'Addition - Recouvrement
    for i = 2:K - 1 %On enlève la première et le dernière trame
        
        %Découpe en trames
        start_idx = (i-1)*ecart + 1;
        stop_idx = start_idx + NFFT - 1;
    
        %Créer une fenêtre de la taille de la trame x
        window = type_fenetre(length(trames_restaurees(:,1)));
        
        %Fenetrage de la fenetre
        trame_windowed = trames_restaurees(:,i) .* window;
        
        %Reconstitution
        %s_rebuilt(start_idx:stop_idx) = trame_windowed; 
        s_rebuilt(start_idx:stop_idx) = s_rebuilt(start_idx:stop_idx) + trame_windowed;
    end
    % Finalisation: Assurer la bonne taille et ajuster (normalisation par la fenêtre)
    signal_rehausse = s_rebuilt(1:L); % Troncature à la longueur originale
    
    % Normalisation simple (souvent nécessaire quand on utilise OLA)
    % Vous pouvez commenter cette ligne si le son est trop fort, ou si vous
    % implémentez une normalisation plus sophistiquée de la fenêtre.
    signal_rehausse = signal_rehausse / (1/overlap);
end