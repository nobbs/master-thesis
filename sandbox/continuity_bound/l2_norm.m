function y = l2_norm( fun, xspan )
    y = sqrt(quadgk(@(x) abs(fun(x)).^2, xspan(1), xspan(2)));
end

