function [coeffs, y] = entw2(f, grid, N, weights, L, variante, to_plot)
% if nargin == 4
%     variante = 0;
%     to_plot = 0;
% end

coeffs = zeros(N, 1);

% Koeffizienten berechnen
if variante == 0
    % Entwicklung in Sinusfunktionen
    for k = 1:N
        fc        = L * weights(k)^2 / 2;
        coeffs(k) = trapez(weights(k) * f .* sin(pi * k * grid / L), grid, 0) / fc;
    end
elseif variante == 1
    % Entwicklung in Cosinusfunktionen
    for k = 1:N
        kk = k - 1;
        if kk == 0
            fc = L * weights(k)^2;
        else
            fc = L * weights(k)^2 / 2;
        end
        coeffs(k) = trapez(weights(k) * f .* cos(pi * kk * grid / L), grid, 0) / fc;
    end
elseif variante == 2
    % Fourierentwicklung
    for k = 1:N
        kk = k - 1;
        if kk == 0
            fc = L * weights(1)^2 / 2;
            coeffs(1) = trapez(f, grid, 0) / fc;
        else
            fc1 = L * weights(2*kk)^2 / 2;
            fc2 = L * weights(2*kk+1)^2 / 2;
            coeffs(2*kk)   = trapez(weights(2*kk) * f .* cos(2 * pi * kk * grid / L), grid, 0)    / fc1;
            coeffs(2*kk+1) = trapez(weights(2*kk +1) * f .* sin(2 * pi * kk * grid / L), grid, 0) / fc2;
        end
    end
elseif variante == 3
    % Entwicklung in Legendre-Polynome
    for k = 1:N
        kk = k - 1;
        fc = L * weights(k)^2 / (2*kk + 1);
        coeffs(k) = trapez(weights(k) * f .* legendre_P(grid / L, kk), grid, 0) / fc;
    end
end

% Funktion wieder rekonstruieren
if variante == 0
    % Entwicklung in Sinusfunktionen
    y = zeros(size(grid, 1), size(grid, 2));
    for k = 1:N
        y = y + coeffs(k) * weights(k) * sin(pi * k * grid / L);
    end
elseif variante == 1
    % Entwicklung in Cosinusfunktionen
    y = zeros(size(grid, 1), size(grid, 2));
    for k = 1:N
        kk = k - 1;
        y  = y + coeffs(k) * weights(k) * cos(pi * kk * grid / L);
    end
elseif variante == 2
    % Fourierentwicklung
    y = zeros(size(grid, 1), size(grid, 2));
    for k = 1:N
        kk = k - 1;
        if kk == 0
            y = y + coeffs(1) / 2;
        else
            y = y + coeffs(2 * kk) * weights(2 * kk) * cos(2 * pi * kk * grid / L) + coeffs(2 * kk + 1) * weights(2 * kk + 1) * sin(2 * pi * kk * grid / L);
        end
    end
elseif variante == 3
    % Entwicklung in Legendre-Polynome
    y = zeros(size(grid, 1), size(grid, 2));
    for k = 1:N
        kk = k - 1;
        y  = y + coeffs(k) * weights(k) * legendre_P(grid / L, kk);
    end
end

if to_plot 
    % Debugging: Plot...
    figure(1);
    semilogy(abs(coeffs), 'o:');
    ylim([1e-12, 1e2]);
    title('Koeffizienten');
    figure(2);
    plot(grid, f, grid, y);
    title('Funktionen');
end

end
