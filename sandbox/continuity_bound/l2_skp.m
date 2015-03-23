function y = l2_skp( fun1, fun2, xspan )
    y = quadgk(@(x) fun1(x) .* fun2(x), xspan(1), xspan(2));
end

