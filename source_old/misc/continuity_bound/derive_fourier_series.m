function [ dcoeffs ] = derive_fourier_series( N, coeffs )
    dcoeffs = zeros(2*N+1, 1);
    dcoeffs(2*(1:N)) = - (1:N)' .* coeffs(2*(1:N)+1);
    dcoeffs(2*(1:N) + 1) = (1:N)' .* coeffs(2*(1:N));
    dcoeffs = 2 * pi * dcoeffs;
end

