function schranke(vals, grid, L, var, N, ep)
% SCHRANKE Berechnet die Koeffizienten einer Funktion bzgl. einer Entwicklung.
% Variante:
% 0: Sinus-Entwicklung
% 1: Cosinus-Entwicklung
% 2: Fourier-Entwicklung
% 3: Legendre-Entwicklung (TODO!)
% Parameter:
% ep  Epsilon
% N  Anzahl Koeffizienten (bei Fourier 2*N + 1)!

if nargin == 3
    var = 0;
    N = 50;
    ep = 0.01;
elseif nargin == 4
    N = 50;
    ep = 0.01;
elseif nargin == 5
    ep = 0.01;
end

if var == 0
    % Schranke bei Sinus-Entw.
    c = pi^(4 + ep) / (2 * L^3 * zeta(2 + ep))
elseif var == 1
    % Schranke bei Cosinus-Entw.
    c = pi^(4 + ep) / (L^3 * (pi^(2 + ep) + zeta(2 + ep)))
elseif var == 2
    % Schranke bei Fourier-Entw.
    c = pi^(4 + ep) / (L^3 * (pi^(2 + ep) + 3 * zeta(2 + ep)))
elseif var == 3
    % Schranke bei Legendre-Entw.
    % TODO!
    c = 1e-16;
end

% Gewichte
if var == 0
    weights = zeros(N, 1);
    for k = 1:N
        weights(k) = 1 / (pi * k)^(1 + ep);
    end
elseif var == 1
    weights = zeros(N, 1);
    for k = 1:N
        kk = k - 1;
        if kk == 0
            weights(1) = 1;
        else
            weights(k) = 1 / (pi * kk)^(1 + ep);
        end
    end
elseif var == 2
    weights = zeros(2 * N + 1, 1);
    for k = 1:N
        kk = k - 1;
        if kk == 0
            weights(1) = 1;
        else
            weights(2*kk)   = 1 / (pi * kk)^(1 + ep);
            weights(2*kk+1) = 1 / (pi * kk)^(1 + ep);
        end
    end
elseif var == 3
    weights = ones(N, 1);
end

% Vorbereitung
grid = [grid, L];
f    = [vals; vals(1)] - vals(1);

% Koeffizienten bestimmen
[coeffs, y] = entwicklung(f.', grid, N, weights, L, var, 0);

figure(1)
plot(grid, f, grid, y);
title('Funktion und zugehörige Entwicklung');
figure(2)
M = length(coeffs);
plot(1:M, abs(coeffs), 'o:', 1:M, c * ones(M, 1));
title('Koeffizienten und Schranke');
end
