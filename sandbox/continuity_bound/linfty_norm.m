function [ y ] = linfty_norm( fun, steps )
    if nargin == 1
        steps = 100;
    end
    
    y = max(abs(fun(linspace(0, 1, steps))));
end

