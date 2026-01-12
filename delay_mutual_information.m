function [tau_opt, I_tau] = delay_mutual_information(x, maxTau, nBins)
% DELAY_MUTUAL_INFORMATION
% Wyznacza optymalne opóźnienie czasowe na podstawie
% średniej informacji wzajemnej (Fraser & Swinney, 1986)
%
% WEJŚCIE:
%   x       - jednowymiarowy szereg czasowy (wektor)
%   maxTau  - maksymalne opóźnienie (w próbkach)
%   nBins   - liczba binów histogramu (np. 32 lub 64)
%
% WYJŚCIE:
%   tau_opt - optymalne opóźnienie (pierwsze minimum AMI)
%   I_tau   - wartości informacji wzajemnej dla tau = 1:maxTau

x = x(:);
N = length(x);

% normalizacja (opcjonalna, ale stabilizuje histogramy)
x = (x - mean(x)) / std(x);

I_tau = zeros(maxTau,1);

for tau = 1:maxTau
    x1 = x(1:N-tau);
    x2 = x(1+tau:N);

    % histogramy brzegowe
    [P1, edges1] = histcounts(x1, nBins, 'Normalization', 'probability');
    [P2, edges2] = histcounts(x2, nBins, 'Normalization', 'probability');

    % histogram łączny
    P12 = histcounts2(x1, x2, edges1, edges2, ...
                      'Normalization', 'probability');

    % usuwanie zer (log(0))
    epsVal = 1e-12;
    P1(P1 == 0)   = epsVal;
    P2(P2 == 0)   = epsVal;
    P12(P12 == 0) = epsVal;

    % informacja wzajemna
    I = 0;
    for i = 1:nBins
        for j = 1:nBins
            I = I + P12(i,j) * log2( P12(i,j) / (P1(i)*P2(j)) );
        end
    end

    I_tau(tau) = I;
end

% --- wybór opóźnienia: pierwsze minimum lokalne ---
tau_opt = NaN;
for k = 2:maxTau-1
    if I_tau(k) < I_tau(k-1) && I_tau(k) < I_tau(k+1)
        tau_opt = k;
        break;
    end
end

end
