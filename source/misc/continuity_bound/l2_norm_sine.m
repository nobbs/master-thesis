function [ y ] = l2_norm_sine( coeffs, L )
if nargin == 1
    L = 1;
end

    y = sqrt((L / 2) * sum(coeffs.^2));
end

