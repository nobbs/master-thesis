function [ y ] = l2_norm_fourier( coeffs )
    y = sqrt(coeffs(1)^2 / 4 + sum(coeffs(2:end).^2) / 2);
end

