function y = fourier_series(x, N, coeffs)
    y = coeffs(1) / 2;
%     for k = 1:N
%         y = y + coeffs(2 * k) * sin(2 * pi * k * x) + coeffs(2 * k + 1) * cos(2 * pi * k * x); 
%     end
    size(cos(2 * pi * (1:N) .* x))
    y = y + sum(coeffs(2*(1:N))' * cos(2 * pi * (1:N) .* x));
    y = y + sum(coeffs(2*(1:N) + 1)' * sin(2*pi* (1:N) .* x));
end