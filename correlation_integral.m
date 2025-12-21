function [r_vals, Cr] = correlation_integral(signal, m, tau, r_vals)
    % Funkcja oblicza całkę korelacyjną C(r) dla zadanego wymiaru m i opóźnienia tau
    % sygnał - wejściowy szereg czasowy
    % m - wymiar zanurzenia
    % tau - opóźnienie czasowe
    % r_vals - wektor promieni r, dla których liczymy C(r)

    N = length(signal);
    % Liczba punktów w zrekonstruowanej przestrzeni
    M = N - (m-1)*tau;
    
    % Rekonstrukcja przestrzeni fazowej (macierz M x m)
    X = zeros(M, m);
    for i = 1:m
        X(:, i) = signal((1:M) + (i-1)*tau);
    end
    
    % Obliczanie macierzy odległości (uproszczone dla wydajności - pdist)
    % Jeśli masz mało RAMu, można użyć podzbioru punktów (np. 1000-2000)
    if M > 2000
        indices = round(linspace(1, M, 2000));
        X_sub = X(indices, :);
        M_eff = length(indices);
    else
        X_sub = X;
        M_eff = M;
    end
    
    dists = pdist(X_sub);
    
    % Obliczanie C(r) dla każdego zadanego progu r
    Cr = zeros(size(r_vals));
    total_pairs = M_eff * (M_eff - 1) / 2;
    
    for k = 1:length(r_vals)
        Cr(k) = sum(dists < r_vals(k)) / total_pairs;
    end
end