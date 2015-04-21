function [B, ufun] = start(num_M, num_Q, u0, t_plot)
% Galerkin-Verfahren für eindimensionale Raum-Zeit-Variationsformulierung.
% Das Variationsproblem $b(u, v) = f(v)$ mit $u \in \mathcal X$ und $v \in \mathcal Y$ wird durch
% Galerkin-Projektion auf endlichdimensionale Unterräume $\mathcal X_N$ und $\mathcal Y_N$
% approximiert.

if nargin == 3
    t_plot = 0;
end

% Daten und Einstellungen laden
[data, opt, omega] = setup(num_M, num_Q, u0);

% Gewichtungsmatrizen zur Normierung der Steifigkeitsmatrix und des Lastvektors bestimmen
[NX, NY] = normierungsmatrizen(opt);

% Massematrizen für die diskreten X- und Y-Normen bestimmen
[MX, MY] = diskrete_normen(opt);

% Steifigkeitsmatrix aus den berechneten Werten erzeugen und normieren
B = steifigkeitsmatrix(data, omega, opt);
B = NY * B * NX;

% Lastvektor bestimmen und normieren
F = lastvektor(data, opt);
F = NY * F;

%% Das lineare Gleichungssystem $B u = F$ lösen und die Lösung aufbauen.
% Dazu wird zunächst das LGS gelöst.
u = B \ F;

% Und anschließend der daraus resultierende Koeffizientenvektor verwendet um die Lösung zu
% rekonstruieren.
ufun = @(t, x) solution(t, x, opt, u);

% Und weil's so toll ist, plotten wir die Lösung bei Bedarf auch noch.
if t_plot
    figure(2)
    tgrid = linspace(data.tspan(1), data.tspan(2), 100);
    xgrid = linspace(data.xspan(1), data.xspan(2), 50);
    [T, X] = meshgrid(tgrid, xgrid);
    mesh(T, X, ufun(T, X));
end

end
