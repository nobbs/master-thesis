function [ funs ] = sinus_functions( N )
    funs = cell(N, 1);
    
    for k = 1:N
        funs{k} = @(x) sin(2 * k * pi * x);
    end
end

