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
[B, uf] = start(30, 30, u0);
figure();

mesh(T, X, abs(uf(T, X) - exf(20)));
title('Fehlerplot')
xlabel('t')
ylabel('x')

%%
mesh(T, X, uf(T, X));
