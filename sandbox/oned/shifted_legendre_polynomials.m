function [ funs ] = shifted_legendre_polynoms( N )
    funs = cell(N + 1, 1);
    
    for n = 0:N
        P = @(x) 0 * x;
        for k = 0:n
            P = @(x) P(x) + nchoosek(n, k) * nchoosek(n + k, k) * (-x).^k;
        end
        funs{n + 1} = @(x) (-1)^n * P(x);
    end
end

