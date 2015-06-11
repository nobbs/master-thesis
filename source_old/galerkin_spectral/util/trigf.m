function y = trigf(x, k)
    if k == 0
        y = ones(size(x, 1), size(x, 2));
    elseif k > 0
        y = sin(2 * pi * k * x);
    else
        y = cos(2 * pi * (-k) * x);
    end
end
