tspan = [0 1];
xspan = [0 1];

num_M = 15;
num_Q = 10;

% Anzahl Basisfunktionen für X
num_X_j = num_M;
num_X_k = num_Q + 1;

% Anzahl Basisfunktionen für Y
num_Y_l = num_M;
num_Y_m = num_Q;
num_Y_n = num_M;

% Dimension von X und Y
dim_X = num_X_j * num_X_k;
dim_Y = num_Y_l * num_Y_m + num_Y_n;

assert(dim_X == dim_Y)

% Legendre-Polynome bestimmen
legendre_polys = shifted_legendre_polynomials(max(num_X_k, num_Y_m) + 1);
legendre_polys_derivative = shifted_legendre_polynomials_derivative(max(num_X_k, num_Y_m) + 1);

% X_h = \span{ \sin(j * \pi * x) * P_k(t): j = 1..num_m, k = 1..(num_q + 1) }
% Y_h = \span{ ( \sin(l * \pi * x) * P_m, \sin(n * pi * x) ) : \nu = 1..num_m, \mu = 1..num_q, j = 1..num_m }

% \omega und Anfangswert
w = @(x) 1 + 0 .* x;
u0 = @(x) sin(pi*x);
g = @(t, x) 0 .* t + 0 .* x;


% Steifigkeitsmatrix aufsetzen
idx = [];
idy = [];
entry = [];
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
                    tmp2 = ((pi * j)^2 / (2 * (k - 1) + 1)) / 2;
                end;

                tmp3 = 0;
                if k == m
                    tmp3 = (1 / (2 * (k - 1) + 1)) * quadgk(@(x) w(x) .* sin(pi*j*x) .* sin(pi*l*x), xspan(1), xspan(2));
                end;
                
                x_pos = (j - 1) * num_X_k + k;
                y_pos = (l - 1) * num_Y_m + m;

                value = tmp1 + tmp2 + tmp3;
                if value ~= 0
                    idx(end + 1) = x_pos;
                    idy(end + 1) = y_pos;
                    entry(end + 1) = value;
                end
            end
        end
        
        for n = 1:num_Y_n
            if j == n
                tmp4 = (-1)^(k - 1) / 2;
                x_pos = (j - 1) * num_X_k + k;
                y_pos = num_Y_l * num_Y_m + n;
                
                idx(end + 1) = x_pos;
                idy(end + 1) = y_pos;
                entry(end + 1) = tmp4;
            end;
        end
    end
end
B = sparse(idy, idx, entry, dim_Y, dim_X);

% Lastvektor
F = sparse(dim_Y, 1);
for l = 1:num_Y_l
    for m = 1:num_Y_m
        pos = (l - 1) * num_Y_m + m;
        F(pos) = integral2(@(t, x) g(t, x) .* sin(pi*l*x) .* legendre_polys{m}(t), tspan(1), tspan(2), xspan(1), xspan(2));
    end
end
for n = 1:num_Y_n
    pos = num_Y_l * num_Y_m + n;
    F(pos) = integral(@(x) u0(x) .* sin(pi*n*x), xspan(1), xspan(2));
end

% LGS lösen
u = B \ F;

% Lösung zusammenbauen
ufun = @(t, x) 0;
for j = 1:num_X_j
    for k = 1:num_X_k
        ufun = @(t, x) ufun(t, x) + u((j - 1) * num_X_k + k) * sin(pi*j*x) .* legendre_polys{k}(t);
    end
end

if t_plot
    figure(2)
    tgrid = linspace(tspan(1), tspan(2), 50);
    xgrid = linspace(xspan(1), xspan(2), 50);
    [T, X] = meshgrid(tgrid, xgrid);
    mesh(T, X, ufun(T, X));
end

