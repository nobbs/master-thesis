function cosine_series(grid, values, N)
%COSINE_SERIES Rekonstruktion einer Funktion mittels Cosinus-Reihe

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

% Ansatzfunktionen
function y = phi(x, kk)
    if kk == 0
        y = ones(size(x, 1), size(x, 2));
    else
        y = ak(kk + 1) * cos(pi * (kk) * x / L);
    end
end

% Koeffizienten bestimmen
sigma = zeros(N, 1);
for k = 1:N
    kk = k - 1;
    if kk == 0
        ct = cumtrapz(grid, values);
        sigma(k) = ct(end) / L;
    else
        ct = cumtrapz(grid, values .* phi(grid, kk));
        sigma(k) = ct(end) * 2 / (ak(k)^2 * L);
    end
end

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
