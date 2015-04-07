function [ y] = h1_norm_sine( coeffs, L )
    if nargin == 1
        L = 1;
    end;
    
    K = length(coeffs);
    dcoeffs = pi * (1:K)' .* coeffs / L;
    
    y1 = (L / 2) * sum(coeffs.^2);
    y2 = (L / 2) * sum(dcoeffs.^2);
    y = sqrt(y1 + y2);
end

