function [ funs ] = shifted_legendre_polynoms_derivative( N )
    funs = cell(N + 1, 1);
    
    for n = 0:N
        P = @(x) 0 * x;
        for k = 1:n
            P = @(x) P(x) + (-1)^k * k * nchoosek(n, k) * nchoosek(n + k, k) * (-x).^(k-1);
        end
        funs{n + 1} = @(x) (-1)^n * P(x);
    end
end

