function DSP = daniell(x, M)
    % DSP calculée avec la méthode de Daniell. M est la taille de fenetre,
    % influence sur la résolution et la variance du signal. 
    N = length(x);
    DSP_init = (abs(fft(x)).^2) / N;
    DSP = zeros(1, N);
    for i=1:N
        i_deb = max(1, i - M);
        i_fin = min(N, i + M);
        DSP(i) = mean(DSP_init(i_deb:i_fin));
    end
end