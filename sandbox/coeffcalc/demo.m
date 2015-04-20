%% Laden
load('omegas');

% Variante:
% 0: Sinus-Entwicklung
% 1: Cosinus-Entwicklung
% 2: Fourier-Entwicklung

%% Auf Omega = [0, 1] skaliert:
schranke(occh.o1{1}, xgrid / data.L, 1, 3);

%% Nicht skaliert, d.h. Omega = [0, L]:
schranke(occh.o1{1}, xgrid, data.L, 0);
