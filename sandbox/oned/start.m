tspan = [0 1];
xspan = [0 1];

% Anzahl Basisfunktionen f�r X
num_X_j = 15;
num_X_k = 15;

% Anzahl Basisfunktionen f�r Y
num_Y_l = 15;
num_Y_m = 15;
num_Y_n = 15;

% Dimension von X und Y
dim_X = num_X_j * num_X_k;
dim_Y = num_Y_l * num_Y_m * num_Y_n;

% Legendre-Polynome bestimmen
legendre_polys = shifted_legendre_polynomials(max(num_X_k, num_Y_m) + 1);
legendre_polys_derivative = shifted_legendre_polynomials_derivative(max(num_X_k, num_Y_m) + 1);

% X_h = \span{ \sin(j * \pi * x) * P_k(t): j = 1..num_m, k = 1..(num_q + 1) }
% Y_h = \span{ ( \sin(l * \pi * x) * P_m, \sin(n * pi * x) ) : \nu = 1..num_m, \mu = 1..num_q, j = 1..num_m }

% Steifigkeitsmatrix aufsetzen
B = sparse(dim_Y, dim_X);

% \omega und Anfangswert
w = @(x) 0 + 0 .* x;
u0 = @(x) x;

for j = 1:num_X_j
    for k = 1:num_X_k
        for l = 1:num_Y_l
            for m = 1:num_Y_m
                tmp1 = 0;
                if j == l
                    tmp1 = quadgk(@(t) legendre_polys_derivative{k}(t) .* ...
                            legendre_polys{m}(t), tspan(1), tspan(2)) / 2;
                end;

                tmp2 = 0;
                if j == l && k == m
                    tmp2 = ((2 * pi * j)^2 / (2 * (k - 1) + 1)) / 2;
                end;

                tmp3 = 0;
                if k == m
                    tmp3 = 1 / (2 * (k - 1) + 1) * quadgk(@(x) w(x) .* sin(2*pi*j*x) .* sin(2*pi*l*x), xspan(1), xspan(2));
                end;

                for n = 1:num_Y_n
                    y_pos = (j - 1) * num_X_k + k;
                    x_pos = (l - 1) * (num_Y_m * num_Y_n) + (m - 1) * num_Y_n + n;

                    tmp4 = 0;
                    if (j == n)
                        tmp4 = (-1)^(k - 1) / 2;
                    end;

                    B(x_pos, y_pos) = tmp1 + tmp2 + tmp3 + tmp4;
                end
            end
        end
    end
end

% Lastvektor
F = sparse(dim_Y, 1);
for l = 1:num_Y_l
    for m = 1:num_Y_m
        for n = 1:num_Y_n
            pos = (l - 1) * (num_Y_m * num_Y_n) + (m - 1) * num_Y_n + n;
            F(pos) = quadgk(@(x) u0(x) .* sin(2*pi*n*x), xspan(1), xspan(2));
        end
    end
end

u = B \ F;

% L�sung zusammenbauen
ufun = @(t, x) 0;
for j = 1:num_X_j
    for k = 1:num_X_k
        ufun = @(t, x) ufun(t, x) + u((j - 1) * num_X_k + k) * sin(2*pi*j*x) * legendre_polys{k}(t);
    end
end

