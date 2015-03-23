function y = generate_sine_series(x, coeffs)
    N = length(coeffs);
    y = zeros(size(x,1), size(x,2));
    for k = 1:N
        y = y + coeffs(k) * sin(pi * k * x); 
    end
end
