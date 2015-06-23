function y = omega_fourier(x, epsilon, N, coeffs)
%OMEGA_SINUS Berechnet die Reihe `\sum_{j=0}^{N} \sigma_j \varphi_j`, wobei
% `\sigma` = coeffs und die Ansatzfunktionen gewählt sind als
% `\varphi_j(x) = \frac{K}{(j\pi)^{1 + \epsilon}} \sin(j \pi x)`.
    N_half = (N - 1) / 2;
    % y = coeffs(1) * ones(size(x, 1), size(x, 2));
    y = zeros(size(x, 1), size(x, 2));
    for k = (-N_half):(N_half)
        if k == 0
            y = y + coeffs(k + N_half + 1);
        else
            y = y + coeffs(k + N_half + 1) * trigf(x, k) / (2 * abs(k) * pi)^(1 + epsilon);
        end
    end
end
