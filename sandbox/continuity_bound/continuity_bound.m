% Berechnet die Stetigkeits-Konstante für $A_l$ mit Ansatzfunktion
% $\phi_l = K * \frac{\sin(\pi * l * x)}{(l * \pi)^{1 + \epsilon})}$.

% Zu testende Ansatzfunktionen
Lspan = 1:25;
% Parameter für Ansatzfunktion
epsilon = 1 / 100;
K       = 5;

% Anzahl Sinuskoeffizienten (Basisfunktionen)
N = 5;

% Anzahl Schleifendurchläufe
M = 100;

% Definitionsbereich
xspan = [0, 1];
values = zeros(M, 1);

% Ergebnisse
bounds = zeros(length(Lspan), 1);
approxs = zeros(length(Lspan), 1);

% Einbettungskonstante
% Cinfty = einbettung();

for idx_l = 1:length(Lspan)
    % Zu testende Ansatzfunktion erzeugen
    l   = Lspan(idx_l);
    phi = @(x) (K / (l * pi)^(1+epsilon)) * sin(l * pi * x);

    % Abschätzung berechenen
    bounds(idx_l) = (4 * K * Cinfty) / (l * pi)^(2 + epsilon);
    % bounds(idx_l) = (4 * K) / (l * pi)^(2 + epsilon);

    % values = zeros(N * M, 1);
    max_val = 0;

    for k = 1:N
        for iter = 1:M
            % gleichverteilte Sinuskoeffizienten generieren und zentrieren
            ucoeffs = 2 * (rand(k, 1) - 0.5);
            vcoeffs = 2 * (rand(k, 1) - 0.5);

            % H^{1}-Normen bestimmen
            unorm = h1_norm_sine(ucoeffs);
            vnorm = h1_norm_sine(vcoeffs);

            % H^{1}-normierte Funktionen aus Sinuskoeffizienten erzeugen
            u  = @(x) generate_sine_series(x, ucoeffs / unorm);
            v  = @(x) generate_sine_series(x, vcoeffs / vnorm);

            % Wert berechnen
            val     = abs(l2_skp(@(x) phi(x) .* u(x), v, xspan));
            max_val = max(val, max_val);
            % values((k - 1) * M + iter) = val;
        end
    end

    % Approximation an Stetigkeitskonstante ermitteln
    approxs(idx_l) = max_val;
end

% Plotten!
figure();
semilogy(Lspan, bounds, Lspan, approxs);
title('Vergleich Abschätzung und experimentelle Stetigkeitskonstante');
xlabel('l, wobei \phi_l(x) = c_l * sin(l \pi x)');
legend('Abschätzung', 'Experimentell');
