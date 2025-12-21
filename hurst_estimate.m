function H = hurst_estimate(signal)
    % Funkcja szacuje wykładnik Hursta metodą R/S (Rescaled Range)
    % signal - wejściowy szereg czasowy
    
    N = length(signal);
    % Dzielimy sygnał na bloki o różnych długościach n
    n_vals = floor(logspace(log10(10), log10(N/4), 10));
    RS = zeros(size(n_vals));
    
    for i = 1:length(n_vals)
        n = n_vals(i);
        m = floor(N/n); % liczba bloków
        rs_blocks = zeros(m, 1);
        
        for j = 1:m
            block = signal((j-1)*n + 1 : j*n);
            mean_b = mean(block);
            % Kumulatywne odchylenie
            Y = cumsum(block - mean_b);
            % Zakres R i odchylenie standardowe S
            R = max(Y) - min(Y);
            S = std(block);
            if S == 0, S = 1; end % Zabezpieczenie przed dzieleniem przez 0
            rs_blocks(j) = R / S;
        end
        RS(i) = mean(rs_blocks);
    end
    
    % Wykładnik Hursta to nachylenie prostej log(n) vs log(RS)
    coeffs = polyfit(log(n_vals), log(RS), 1);
    H = coeffs(1);
    
    % Opcjonalna wizualizacja dopasowania (do sprawozdania)
    %plot(log(n_vals), log(RS), 'o'); hold on; plot(log(n_vals), polyval(coeffs, log(n_vals)));
end