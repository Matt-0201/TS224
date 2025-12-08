function [s_bruit, bbgc_ajusted] = bruitageSignal(s, RSB)
    
    N = length(s);
    bbgc = randn(1,N);                   % Génération d'un bbgc

    Ps = (1/N)*sum(s.^2);                % Puissance du signal
    Pb= (1/N)*sum(bbgc.^2);              % Puissance du bruit

    alpha = sqrt(Ps/(Pb * 10^(RSB/10))); %Ajustement

    bbgc_ajusted = bbgc * alpha;

    s_bruit = s + bbgc_ajusted';
end