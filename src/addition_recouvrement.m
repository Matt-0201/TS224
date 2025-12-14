function s_rebuilt = addition_recouvrement(s, NFFT, type_fenetre, overlap)
    % Fonction qui reconstruit le signal à partir de trames en utilisant la méthode d'addition–recouvrement.
    
    % Entrées :
    % s            : signal d'entrée à découper et reconstruire (vecteur colonne)
    % NFFT         : taille des trames (nombre d'échantillons par trame)
    % type_fenetre : type de fenêtre appliquée aux trames (ex : @hamming, @hann, @blackman, etc)
    % overlap      : taux de recouvrement entre deux trames (entre 0 et 1)
    
    % Sortie :
    % s_rebuilt    : signal reconstruit après fenêtrage et addition–recouvrement    
    
    %% Paramètres
    L = length(s);
    recouvrement = floor(overlap * NFFT); %Nombre de points qui recouvrent morceau du signal
    ecart = NFFT - recouvrement; %Ecart entre chaque morceau
    K = floor((L - NFFT) / ecart) + 1; %nombre total de trames
    s_rebuilt = zeros(L, 1); %signal à reconstruire

    fprintf('Nb trames =  %d : NFFT = %d, ecart = %d\n', K, NFFT, ecart);
    
    %% Boucle d'Addition - Recouvrement
    for i = 2:K - 1 %On enlève la première et le dernière trame
        
        %Découpe en trames
        start_idx = (i-1)*ecart + 1;
        stop_idx = start_idx + NFFT - 1;
        
        %Extraction + zero padding si nécessaire
        if stop_idx <= L
            trame = s(start_idx:stop_idx);
        else %Cas si la dernière trame est incomplete
            trame = zeros(NFFT,1);
            len_restante = L - start_idx + 1;
            trame(1:len_restante) = s(start_idx:L); %Le reste de la trame est remplie de 0 (élement neutre)
        end

        %Créer une fenêtre de la taille de la trame x
        window = type_fenetre(length(trame));
        
        %Fenetrage de la trame
        trame_windowed = trame .* window;
         
        %Reconstitution
        %s_rebuilt(start_idx:stop_idx) = trame_windowed; 
        s_rebuilt(start_idx:stop_idx) = s_rebuilt(start_idx:stop_idx) + trame_windowed;
    end
end

