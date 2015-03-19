function [ funs ] = shifted_legendre_polynoms_derivative( N )
    funs = cell(N + 1, 1);
    funs{1} = @(x) 0 * x;
    funs{2} = @(x) 2 + 0 * x;
    
    for k = 3:N
        n = k - 1;
        funs{k} = @(x) (2*n - 1) / (n - 1) * (2*x - 1) .* funs{k - 1}(x) - (n) / (n - 1) * funs{k-2}(x);
    end
end

