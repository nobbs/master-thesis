function [ y ] = l2_norm_sine( coeffs )
    y = sqrt((1 / 2) * sum(coeffs.^2));
end

