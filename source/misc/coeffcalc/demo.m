%% Laden
load('omegas');

%% Alle omega's plotten
figure();
for k = 1:length(occh.o1)
    plot(xgrid, occh.o1{k});
    drawnow;
    waitforbuttonpress;
end

%% Variante:
% 0: Sinus-Entwicklung
% 1: Cosinus-Entwicklung
% 2: Fourier-Entwicklung

%% Auf Omega = [0, 1] skaliert:
% schranke(occh.o1{end}, xgrid / data.L, 1, 0)

% Nicht skaliert, d.h. Omega = [0, L]:
schranke(occh.o1{4}, xgrid, data.L, 0);
