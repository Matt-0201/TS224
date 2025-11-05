clear 
close all
clc

N = 1024;
sig2 = 1;

bruit = randn(1, N) * sqrt(sig2);

freq = -1:1/N:1;

auto_corr = xcorr(bruit, N);
xcorr_biased = xcorr(bruit, N, 'biased');
xcorr_unbiased = xcorr(bruit, N, 'unbiased');
xcorr_normalized = xcorr(bruit, N, 'normalized');

DSP = abs(fftshift(fft(bruit,N))).^2/N;
f = -1/2:1/N:1/2-1/N;

figure;
hold on;
plot(freq, auto_corr);
plot(freq,xcorr_biased, 'r', LineWidth=4);
plot(freq,xcorr_unbiased, 'g', LineWidth=4);
plot(freq,xcorr_normalized, 'b');
hold off;
legend('Théorique','Biaisé', 'Non biaisé', 'Normalized');
title("Tracé de différents estimateurs de la fonction d'autocorrélation d'un bruit AWGN");

figure;
plot(f, DSP);
title("Densité spectrale de puissance d'un bruit AWGN");

