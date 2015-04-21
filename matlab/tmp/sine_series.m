function sine_series(grid, values, N)
%SINE_SERIES Rekonstruktion einer Funktion mittels Sinus-Reihe

grid = [grid; 5];
values = [values; values(1)];

% Parameter
L = 5;
epsilon = 0.01;

% Koeffizienten der Ansatzfunktionen
ak = zeros(N);
ak(1) = 1;
for k = 2:N
    kk = k - 1;
    ak(k) = 1 / (pi * (kk))^(1 + epsilon);
end

% Matrix aufbauen
A = sparse(N, N);
for k = 1:N
    kk = k - 1;
    % Erste Zeile und erste Spalte
    if mod(kk, 2) == 1
        A(1, k) = 2 * ak(k) * L / (pi * kk);
        A(k, 1) = 2 * ak(k) * L / (pi * kk);
    end
    
    % Hauptdiagonale
    if k == 1
        A(k, k) = L;
    else
        A(k, k) = ak(k)^2 * L / 2;
    end
end

% Ansatzfunktionen
function y = phi(x, kk)
    if kk == 0
        y = ones(size(x, 1), size(x, 2));
    else
        y = ak(kk + 1) * sin(pi * (kk) * x / L);
    end
end

% Vektor
f = zeros(N, 1);
for k = 1:N
    ct = cumtrapz(grid, (values .* phi(grid, k - 1)).');
    f(k) = ct(end);
end

% Koeffizienten bestimmen
sigma = A \ f;

% Reihe auswerten
function y = calc(x)
    y = zeros(size(x, 1), size(x, 2));
    for k = 1:N
        y = y + sigma(k) * phi(x, k - 1);
    end
end

approx = calc(grid);
err = norm(approx - values, 2)
plot(grid, approx, grid, values);
end
