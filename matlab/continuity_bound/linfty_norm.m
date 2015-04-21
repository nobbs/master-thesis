function [ y ] = linfty_norm( fun, steps, L )
    if nargin == 1
        steps = 100;
        L = 1;
    elseif nargin == 2
        L = 1;
    end
    
    y = max(abs(fun(linspace(0, L, steps))));
end

