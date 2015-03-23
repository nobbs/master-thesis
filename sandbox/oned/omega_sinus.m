function y = omega_sinus(x, epsilon, N, coeffs)
%OMEGA_SINUS Berechnet die Reihe $\sum_{j=0}^{N} \sigma_j \varphi_j$, wobei
% \sigma = coeffs und die Ansatzfunktionen gewählt sind als
% $\varphi_j(x) = \frac{K}{(j\pi)^{1 + \epsilon} \sin{j \pi x)$.
    y = coeffs(1) * ones(size(x, 1), size(x, 2));
    for j = 1:N
       y = y + coeffs(j + 1) * sin(pi * j * x) / (j * pi)^(1 + epsilon);
    end
end
