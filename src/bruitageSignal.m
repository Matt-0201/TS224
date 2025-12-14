function [s_bruit, bbgc_ajusted] = bruitageSignal(s, RSB)
    % Fonction qui permet de bruiter un signal de parole avec un BBGC afin d'obtenir un RSB donné.
    
    % Entrées
    % s   : signal de parole original 
    % RSB : rapport signal à bruit souhaité [dB]
    
    % Sorties :
    % s_bruit        : signal de parole bruité avec le RSB imposé
    % bbgc_ajusted   : bruit blanc gaussien centré ajusté
    
    N = length(s);
    bbgc = randn(1,N);                   % Génération d'un bbgc

    Ps = (1/N)*sum(s.^2);                % Puissance du signal
    Pb= (1/N)*sum(bbgc.^2);              % Puissance du bruit

    alpha = sqrt(Ps/(Pb * 10^(RSB/10))); %Ajustement

    bbgc_ajusted = bbgc * alpha;

    s_bruit = s + bbgc_ajusted';
end