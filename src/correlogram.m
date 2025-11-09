function DSP = correlogram(x, M)
    % Estimation de la DSP par la méthode du corrélogramme. 
    % M est la taille de la fenetre utilisé pour calculer la fft sur
    % l'autocorrélation de x.

    x_corr = xcorr(x, "unbiased");
    centre = length(x);

    window = hamming(2*M+1).';

    x_corr_window = x_corr(centre-M:centre+M).*window;
    
    DSP = abs(fftshift(fft(x_corr_window)));
end