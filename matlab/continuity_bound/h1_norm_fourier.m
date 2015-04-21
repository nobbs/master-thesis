function [ y] = h1_norm_fourier( coeffs, dcoeffs )
    y1 = coeffs(1)^2 / 4 + sum(coeffs(2:end).^2) / 2;
    y2 = dcoeffs(1)^2 / 4 + sum(dcoeffs(2:end).^2) / 2;
    y = sqrt(y1 + y2);
end

