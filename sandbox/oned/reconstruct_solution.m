function [ y ] = reconstruct_solution( t, x, num_X_j, num_X_k, u )
%RECONSTRUCT_SOLUTION Baut aus dem berechneten Lösungsvektor eine Lösung.
% Dazu wird mittels der berechneten Koeffizienten eine Linearkombination der
% Ansatzfunktionen gebildet.

y = zeros(size(t, 1), size(t, 2));

for j = 1:num_X_j
    for k = 1:num_X_k
        y = y + u((j - 1) * num_X_k + k) * sin(pi*j*x) .* legendre_P(t, k - 1);
    end
end

end

