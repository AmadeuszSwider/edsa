function Y = create_Y(s, n, T, d)
    % s - sygnał (wektor)
    % n - liczba wierszy macierzy
    % T - przesunięcie (opóźnienie)
    % d - liczba kolumn

    Y = zeros(n, d);
    idx = (0:d-1)*T;
    Y = s((1:n)' + idx);

end
