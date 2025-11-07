function DSP = welch(x, NFFT, Fe)
    %y = Mon_Welch_(x, NFFT, Fe);
              
    nb_points = length(x);
    nb_seg = floor(nb_points/NFFT);
    disp(nb_seg);
    DSP_sum = zeros(1,NFFT);
    for i=1:nb_seg
        seg = x((i-1)*NFFT+1:i*NFFT);
        FFT = fftshift(fft(seg, NFFT));
        DSPi= abs(FFT).^2;
        DSP_sum = DSP_sum + DSPi;
    end

    DSP = (DSP_sum/nb_seg)/NFFT*Fe;
end