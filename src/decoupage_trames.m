function signal_decoupe = decoupage_trames(signal_entree, NFFT, overlap)
    %Fonction qui découpe un signal en K trames
    
    s = signal_entree;
    L = length(s);
    recouvrement = floor(overlap * NFFT); %Nombre de points qui recouvrent morceaux du signal
    ecart = NFFT - recouvrement; %Ecart entre chaque morceau
    K = floor((L - NFFT) / ecart) + 1; %nombre total de trames
    
    %Signal de sortie
    signal_decoupe = zeros(NFFT, K);
    
    for i = 1:K %Pour l'instant, on garde la première et la dernière trame
    
        %Découpe en trames
        start_idx = (i-1)*ecart + 1;
        stop_idx = start_idx + NFFT - 1;
    
        %Extraction + zero padding si nécessaire
        if stop_idx <= L
            signal_decoupe(:, i) = s(start_idx:stop_idx);
        else %Cas si la dernière trame est incomplete
            trame = zeros(NFFT, 1);
            len_restante = L - start_idx + 1;
            trame(1:len_restante) = s(start_idx:L);
            signal_decoupe(:, i) = trame;
        end
    end
end