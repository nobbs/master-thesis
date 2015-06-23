%% Vorbereitung
tspan = [0, 1];
xspan = [0, 1];

tgrid = linspace(tspan(1), tspan(2), 100);
xgrid = linspace(xspan(1), xspan(2), 100);

[T, X] = meshgrid(tgrid, xgrid);

% Exakte Lösung bestimmen.
% k gibt die Anzahl der verwendeten Koeffizienten an.
exf = @(k) exfun(u0, T, X, k);

%% Spektral-Galerkin
[data, opt, omega] = setup(10, 10, u0);
[B, uf] = start(data, opt, omega);
figure();

mesh(T, X, abs(uf(T, X) - exf(20)));
title('Fehlerplot')
xlabel('t')
ylabel('x')

%%
figure()
mesh(T, X, uf(T, X));

%%
figure()
% mesh(T, X, exf(1));
Y = exf(10);
Z = uf(T, X);
plot(Z(:, 1))
