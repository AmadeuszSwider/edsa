%% Projekt 2: Analiza systemów chaotycznych
clear all; close all; clc;

% Parametry globalne
N_target = 10000; % Docelowa liczba próbek 
h_chen = 0.005;

%% 1. Sygnał losowy
s_rand = rand(N_target, 1);

%% 2. Sygnał periodyczny
n_idx = (0:N_target-1)';
s_per = sin(n_idx/10) + cos(n_idx/10);

%% 3. System Lorenza (Zakładając dostępność funkcji FOLorenz)
% Aby mieć 10000 próbek przy kroku 0.01, potrzebujemy TSim = 100
[t_lor, y_lor] = FOLorenz([10 28 8/3], [1 1 1], N_target*0.01, [0.1 0.1 0.1]);
s_lorenz = y_lor(1:N_target, 1);

%% 4. System Chena (całkowity)
% Poprawka: TSim musi być nieco większy, aby wygenerować n+1 punktów
TSim_Chen = N_target * h_chen; 
[t_chen_int, y_chen_int] = FOChen([35 3 28 -7], [1 1 1], TSim_Chen, [-9 -5 14]);

% Pobieramy tyle próbek, ile faktycznie wygenerował system, ale nie więcej niż N_target
n_actual = size(y_chen_int, 1);
s_chen_int = y_chen_int(1:min(n_actual, N_target), 1);

%% 5. System Chena (ułamkowy q=0.95)
[t_chen_frac1, y_chen_frac1] = FOChen([35 3 28 -7], [0.95 0.95 0.95], TSim_Chen, [-9 -5 14]);
[t_chen_frac2, y_chen_frac2] = FOChen([35 3 28 -7], [0.9 0.9 0.9], TSim_Chen, [-9 -5 14]);

s_chen_f2 = y_chen_frac2(1:min(size(y_chen_frac2,1), N_target), 1);
s_chen_f1 = y_chen_frac1(1:min(size(y_chen_frac1,1), N_target), 1);

%% Wizualizacja dla sprawdzenia
figure;
subplot(2,1,1); plot(s_chen_int); title('Chen Całkowity');
subplot(2,1,2); plot(s_chen_f1); title('Chen Ułamkowy q=0.95');

%% 6. Wizualizacja Atraktorów 3D
% Tworzymy figurę z trzema podwykresami dla systemów chaotycznych

figure('Name', 'Atraktory Systemów Chaotycznych', 'NumberTitle', 'off');

% --- Atraktor Lorenza ---
subplot(1, 3, 1);
plot3(y_lor(:,1), y_lor(:,2), y_lor(:,3), 'r');
grid on; axis tight;
title('System Lorenza (Klasyczny)');
xlabel('x(t)'); ylabel('y(t)'); zlabel('z(t)');
view(45, 30); % Ustawienie kąta patrzenia dla lepszej widoczności

% --- Atraktor Chena (Całkowity q=1) ---
subplot(1, 3, 2);
plot3(y_chen_int(:,1), y_chen_int(:,2), y_chen_int(:,3), 'b');
grid on; axis tight;
title('System Chena (q = 1.0)');
xlabel('x(t)'); ylabel('y(t)'); zlabel('z(t)');
view(45, 30);

% --- Atraktor Chena (Ułamkowy q=0.95) ---
subplot(1, 3, 3);
plot3(y_chen_frac1(:,1), y_chen_frac1(:,2), y_chen_frac1(:,3), 'g');
grid on; axis tight;
title('System Chena (q = 0.95)');
xlabel('x(t)'); ylabel('y(t)'); zlabel('z(t)');
view(45, 30);

% Poprawa estetyki wykresów
set(findall(gcf,'-property','FontSize'),'FontSize',10);

%% 7. Estymacja opóźnienia czasowego (Tau) - Autokorelacja
max_lag = 500;
[acf_lor, lags] = autocorr(s_lorenz, 'NumLags', max_lag);

% Szukamy pierwszego przejścia poniżej progu 1/e
threshold = 1/exp(1);
tau_idx = find(acf_lor < threshold, 1);
tau = lags(tau_idx);

figure('Name', 'Analiza Opóźnienia Czasowego');
plot(lags, acf_lor, 'LineWidth', 1.5);
hold on;
line([0 max_lag], [threshold threshold], 'Color', 'r', 'LineStyle', '--');
plot(tau, acf_lor(tau_idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
title(['Autokorelacja dla systemu Lorenza, \tau = ' num2str(tau)]);
xlabel('Przesunięcie (lag)'); ylabel('Wartość ACF');
grid on;

%% 7. Estymacja opóźnienia czasowego (Tau) dla wszystkich sygnałów
signals = {s_rand, s_per, s_lorenz, s_chen_int, s_chen_f1, s_chen_f2};
sig_names = {'Losowy', 'Periodyczny', 'Lorenz', 'Chen (q=1)', 'Chen (q=0.95)', 'Chen (q=0.9)'};
tau_values = zeros(1, length(signals));
max_lag = 500;
threshold = 1/exp(1);

figure('Name', 'Analiza Autokorelacji - Wybór Tau');
for i = 1:length(signals)
    [acf, lags] = autocorr(signals{i}, 'NumLags', max_lag);
    
    % Szukamy pierwszego spadku poniżej 1/e (metoda zalecana)
    idx = find(acf < threshold, 1);
    if isempty(idx)
        tau_values(i) = 1; % Zabezpieczenie dla sygnałów bardzo wolnozmiennych
    else
        tau_values(i) = lags(idx);
    end
    
    subplot(2, 3, i);
    plot(lags, acf); hold on;
    line([0 max_lag], [threshold threshold], 'Color', 'r', 'LineStyle', '--');
    plot(tau_values(i), acf(lags == tau_values(i)), 'ro');
    title([sig_names{i}, ' (\tau=', num2str(tau_values(i)), ')']);
    grid on;
end

%% 8. Analiza nasycenia - Całka korelacyjna i Wymiar zanurzenia de
m_range = 1:5;
r_vals = logspace(-2, 1, 20);
de_values = zeros(1, length(signals)); 
all_D2 = zeros(length(m_range), length(signals)); % Przechowywanie D2 dla tabeli

figure('Name', 'Analiza nasycenia - Wszystkie Sygnały', 'NumberTitle', 'off');

for s = 1:length(signals)
    sig = signals{s};
    tau_s = tau_values(s);
    D2_temp = zeros(length(m_range), 1);
    
    subplot(2, 3, s); hold on;
    for i = 1:length(m_range)
        m = m_range(i);
        [r, Cr] = correlation_integral(sig, m, tau_s, r_vals);
        
        log_r = log(r);
        log_Cr = log(Cr);
        valid = isfinite(log_Cr);
        
        plot(log_r(valid), log_Cr(valid), '-');
        
        % Estymacja nachylenia D2 (wymiar korelacyjny)
        linear_idx = find(log_r > -1.5 & log_r < 0.5 & valid);
        if length(linear_idx) > 1
            p = polyfit(log_r(linear_idx), log_Cr(linear_idx), 1);
            D2_temp(i) = p(1);
        end
    end
    all_D2(:, s) = D2_temp;
    title(sig_names{s});
    xlabel('ln(r)'); ylabel('ln(C(r))');
    grid on;
    
    % Wyznaczanie de (nasycenie) - zgodnie z wytycznymi [cite: 20]
    diffs = diff(D2_temp);
    idx_de = find(diffs < 0.15, 1); % Próg nasycenia
    if isempty(idx_de)
        de_values(s) = max(m_range); 
    else
        de_values(s) = m_range(idx_de);
    end
end

% Wykres nasycenia D2(m) dla wszystkich sygnałów (naprawiony błąd)
figure('Name', 'Krzywe Nasycenia D2(m)');
plot(m_range, all_D2, 's-', 'LineWidth', 1.5);
xlabel('Wymiar zanurzenia m'); ylabel('Wymiar korelacyjny D_2');
legend(sig_names, 'Location', 'best');
title('Porównanie nasycenia wymiaru korelacyjnego');
grid on;


%% 9. Analiza wykładnika Hursta i Tabela Zbiorcza
hurst_values = zeros(1, length(signals));

fprintf('\nObliczanie wykładnika Hursta dla wszystkich sygnałów...\n');
for i = 1:length(signals)
    hurst_values(i) = hurst_estimate(signals{i});
end

% Generowanie tabeli (Wymóg: sekcja 2.2, pkt 50) [cite: 50]
TabelaWynikow = table(sig_names', tau_values', de_values', hurst_values', ...
    'VariableNames', {'Sygnal', 'Tau_Opoznienie', 'de_Zanurzenie', 'Hurst_H'});

disp('--- PODSUMOWANIE PARAMETRÓW REKONSTRUKCJI ---');
disp(TabelaWynikow);

%% 10. Analiza wykładnika Lapunowa (LLE)
% Definicja fs na podstawie kroku h użytym w symulacji
h_lorenz = 0.01; 
fs_lorenz = 1/h_lorenz;

h_chen = 0.005;
fs_chen = 1/h_chen;

% Estymacja LLE (Wymóg projektu [cite: 23])
% Funkcja lyapunovExponent wymaga Toolboxa 'Predictive Maintenance'
try
    lle_lorenz = lyapunovExponent(s_lorenz, fs_lorenz); 
    lle_chen_f = lyapunovExponent(s_chen_f1, fs_chen);
    
    fprintf('\nNajwiększy wykładnik Lapunowa (LLE):\n');
    fprintf('Lorenz: %.4f\n', lle_lorenz);
    fprintf('Chen ułamkowy (q=0.95): %.4f\n', lle_chen_f);
catch
    warning('Funkcja lyapunovExponent nie jest dostępna. Upewnij się, że masz odpowiedni Toolbox.');
end

%% 11. Analiza bifurkacyjna względem rzędu ułamkowego (q)
q_vals = 0.85:0.005:1.0; 
bif_data = [];

fprintf('\nGenerowanie diagramu bifurkacyjnego (wymóg na ocenę > 3.0)...\n');
for q = q_vals
    [~, y_bif] = FOChen([35 3 28 -7], [q q q], 40, [-9 -5 14]);
    % Pomijamy stany przejściowe (pierwsze 60% próbek)
    start_idx = round(size(y_bif,1) * 0.6);
    data_tail = y_bif(start_idx:end, 1); 
    
    % Wyznaczanie maksimów lokalnych
    [pks, ~] = findpeaks(data_tail);
    if ~isempty(pks)
        bif_data = [bif_data; repmat(q, length(pks), 1), pks];
    end
end

figure('Name', 'Diagram Bifurkacyjny');
plot(bif_data(:,1), bif_data(:,2), 'r.', 'MarkerSize', 1);
xlabel('Rząd ułamkowy q'); ylabel('Maksima lokalne x(t)');
title('Przejście do chaosu w ułamkowym systemie Chena');
grid on;

%% Wyświetlenie tabeli parametrów (do sprawozdania)
TabelaWynikow = table(sig_names', tau_values', de_values', hurst_values',...
    'VariableNames', {'Sygnal', 'Tau_Opóźnienie', 'de_Zanurzenie', 'Hurst_H'});
disp(TabelaWynikow);