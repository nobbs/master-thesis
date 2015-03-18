function [ funs ] = sinus_functions_derivative( N )
    funs = cell(N, 1);
    
    for k = 1:N
        funs{k} = @(x) 2 * k * pi * cos(2 * k * pi * x);
    end
end

