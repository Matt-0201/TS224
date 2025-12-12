function signal_rehausse = calcul_fft_trame(signal_en_trames, NFFT, type_fenetre)
    % Fonction réalisant la méthode de la soustraction spectrale
    
    %Paramètres
    K = length(signal_en_trames(1,:));
    signal_rehausse = zeros(NFFT, K);
    
    %Traitement trame à trame
    for i = 1:K %Pour l'instant, on garde la première et la dernière trame
        
        %Transpformée de Fourrier à court terme
        trame = signal_en_trames(:, i); % Extraire la trame actuelle
        
        %Fenêtrage
        trame_windowed = trame .* type_fenetre(NFFT); % Applique la fenêtre à la trame
        fft_trame = fft(trame_windowed); % Calcule la FFT de la trame fenêtrée
        signal_rehausse(:,i) = abs(fft_trame);
    end

    % window = type_fenetre(NFFT);
    % trame_windowed = trame_1 .* window;
    % fft_trame_1 = fft(trame_windowed);
end