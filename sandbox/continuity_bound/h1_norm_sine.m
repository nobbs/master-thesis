function [ y] = h1_norm_sine( coeffs )
    K = length(coeffs);
    dcoeffs = pi * (1:K)' .* coeffs;
    
    y1 = (1 / 2) * sum(coeffs.^2);
    y2 = (1 / 2) * sum(dcoeffs.^2);
    y = sqrt(y1 + y2);
end

