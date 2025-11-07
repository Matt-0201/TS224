function DSP = bartlett(x, NFFT)
% DSP calculée avec la méthode de Barlett
    N = length(x);
    N_zeros = 2^nextpow2(NFFT);
    x = [x zeros(1, N_zeros - NFFT)]; % 0 padding si NFFT n'est pas une puissance de 2
    K = N_zeros/NFFT;
    disp(K);
    DSP = zeros(1,N_zeros);
    M = N_zeros / K;
    disp(M);
    
    for i=1:K
        seg = x((i-1)*M+1:i*M);
        DSP_seg = (abs(fftshift(fft(seg,NFFT))).^2)/M;
        DSP((i-1)*M+1:i*M) = mean(DSP_seg((i-1)*M+1:i*M));
        %DSP = DSP + mean(DSP_seg);
    end
end

