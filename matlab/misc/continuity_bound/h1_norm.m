function y = h1_norm( fun, dfun, xspan )
    y1 = quadgk(@(x) abs(fun(x)).^2, xspan(1), xspan(2));
    y2 = quadgk(@(x) abs(dfun(x)).^2, xspan(1), xspan(2));
    y = sqrt(y1 + y2);
end

