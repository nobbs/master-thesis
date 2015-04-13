function [B] = steifigkeitsmatrix(data, omega, opt)
% STEIFIGKEITSMATRIX Berechnet die Steifigkeitsmatrix.

% Listen für die Erstellung der Sparse-Matrix $B$ erstellen.
idx = zeros(opt.dim.X, 1);
idy = zeros(opt.dim.X, 1);
val = zeros(opt.dim.X, 1);
ctr = 1;

% Über alle Basisfunktionen iterieren für $\mathcal X_N$ und $\mathcal Y_N$
% iterieren. Dabei nutzen wir die Kreuzprodukt-Struktur von $\mathcal Y_N$ aus,
% indem wir die beiden Komponenten einzeln betrachten.
for j = 1:opt.num.Xj
    for k = 1:opt.num.Xk
        % Nur die erste Komponente von $v$, das heißt $v = (v_1, 0)$.
        for l = 1:opt.num.Yl
            for m = 1:opt.num.Ym
                % Die Bilinearform $b(\cdot, \cdot)$ wurde in die einzelnen
                % Terme zerlegt, so dass möglichst viel bereits per Hand
                % vorberechnet werden kann.

                % $\int_T \int_\Omega u_t(t) v_1(t) dx dt$
                term1 = 0;
                if opt.variante_term1 == 0
                    % Variante 1: Integral mittels numerischer Quadratur
                    % approximieren.
                    if j == l
                        term1 = integral(@(t) legendre_dP(t, k - 1) .* ...
                                           legendre_P(t, m - 1), ...
                                         tspan(1), tspan(2)) / 2;
                    end;
                elseif opt.variante_term1 == 1
                    % Variante 2: Ausnutzen, dass die Legendre-Polynome eine
                    % Orthogonalbasis für die Polynome darstellen. Damit
                    % vereinfacht sich alles zu...
                    if j == l
                        if k == 2 && m == 1
                            term1 = 1;
                        elseif k > m && mod(k, 2) ~= mod(m, 2)
                            term1 = 1 / 2;
                        end
                    end
                else
                    assert(false);
                end

                % $\int_T \int_\Omega data.c_D \nabla u(t) \nabla v_1(t) dx dt$
                term2 = 0;
                if j == l && k == m
                    term2 = data.c_D * ((pi * j)^2 / (2 * (k - 1) + 1)) / 2;
                end;

                % $\int_T \int_\Omega \omega u(t) v_1(t) dx dt$
                if opt.variante_term3 == 0
                    % Variante 1: Integral mittels numerischer Quadratur
                    % approximieren.
                    term3 = 0;
                    if k == m
                        term3 = (1 / (2 * (k - 1) + 1)) * ...
                            integral(@(x) omega(x) .* sin(pi*j*x) .* ...
                                       sin(pi*l*x), xspan(1), xspan(2));
                    end;
                elseif opt.variante_term3 == 1
                    % Variante 2: Integral vorher explizit lösen und hier nur
                    % noch auswerten. Problem: das Integral vom Produkt dreier
                    % Sinusfunktionen ist alles andere als schön. Dadurch
                    % ergeben sich hier einige Fallunterscheidungen...
                    term3 = 0;
                    if k == m
                        % Konstanter Term ist noch einfach
                        if j == l
                            term3 = omega.sigmas(1) / 2;
                        end
                        % Die Sinusfunktionen dagegen nicht mehr. Hier wurde
                        % stark vereinfacht und alle Fälle, die automatisch
                        % Null ergeben direkt weggelassen
                        for iter = 1:omega.N
                            tmp = 0;
                            if mod(iter + j + l, 2) == 1
                                denom = (iter - j - l) * (iter + j - l) * ...
                                        (iter - j + l) * (iter + j + l) * pi;
                                if denom ~= 0
                                    tmp = (- 4 * iter * j * l) / denom;
                                end
                            end
                            term3 = term3 + (omega.sigmas(iter + 1) / ...
                                             (iter * pi)^(1 + data.eps)) * tmp;
                        end
                    end;
                    term3 = (1 / (2 * (k - 1) + 1)) * term3;
                else
                    assert(false);
                end;

                % Der berechnete Wert wird natürlich nur dann abgespeichert,
                % wenn er ungleich Null ist.
                value = term1 + term2 + term3;
                if value ~= 0
                    x_pos = (j - 1) * opt.num.Xk + k;
                    y_pos = (l - 1) * opt.num.Ym + m;

                    idx(ctr) = x_pos;
                    idy(ctr) = y_pos;
                    val(ctr) = value;
                    ctr = ctr + 1;
                end
            end
        end

        % Und jetzt die zweite Komponente von $v$, das heißt $v = (0, v_2)$.
        for n = 1:opt.num.Yn
            if j == n
                tmp4 = (-1)^(k - 1) / 2;

                x_pos = (j - 1) * opt.num.Xk + k;
                y_pos = opt.num.Yl * opt.num.Ym + n;

                idx(ctr) = x_pos;
                idy(ctr) = y_pos;
                val(ctr) = tmp4;
                ctr = ctr + 1;
            end;
        end
    end
end

B = sparse(idy, idx, val, opt.dim.Y, opt.dim.X);
end
