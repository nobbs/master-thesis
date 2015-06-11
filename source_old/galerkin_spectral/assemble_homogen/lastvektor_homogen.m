function [F] = lastvektor_homogen(data, opt)
% LASTVEKTOR Berechnet den Lastvektor zu geg. Anfangs- und Randbedingungen.

f = zeros(opt.dim.Y, 1);

% Wie zuvor betrachten wir die beiden Komponenten von $\mathcal Y_N$ getrennt.
% Zunächst also $v = (v_1, 0)$.
for ldx = 1:opt.num.Yl
    for mdx = 1:opt.num.Ym
      mm = mdx - 1;
      L = data.xspan(2);
      val = integral2(@(t, x) data.g(t, x) .* sin(pi*ldx*x / L) .* ...
                          legendre_P_shifted(t, mm, data.tspan), ...
                        data.tspan(1), data.tspan(2), ...
                        data.xspan(1), data.xspan(2));
      pos    = (ldx - 1) * opt.num.Ym + mdx;
      f(pos) = val;
    end
end

% Und nun $v = (0, v_2)$.
for ndx = 1:opt.num.Yn
    val = integral(@(x) data.u0(x) .* sin(pi*ndx*x / L), ...
                     data.xspan(1), data.xspan(2));
    pos    = opt.num.Yl * opt.num.Ym + ndx;
    f(pos) = val;
end

F = sparse(f);

end
