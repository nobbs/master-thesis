%% FEM-Galerkin
start;

%% Exakte Lösung bestimmen.
% k gibt die Anzahl der verwendeten Koeffizienten an.
exf = @(k) exfun(u0, T, X, k);

%%
figure();
ep = 1e-4;
er = abs(upad - exf(40));
eridx = find(er > ep);
er(eridx) = 0;

mesh(T, X, er);
title('Fehlerplot')
xlabel('t')
ylabel('x')

%%
mesh(T, X, upad)

%%
plot(upad(2, :))
