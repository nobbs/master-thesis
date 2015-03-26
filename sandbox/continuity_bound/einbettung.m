function [ y ] = einbettung()
% Berechnet die Einbettungskonstante für die Einbettung $H^{1}_{0}(\Omega)
% \subset L_{\infty}(\Omega)$ anhand von Sinusfunktionen mit zufällig
% gleichverteilten Koeffizienten.

% Anzahl Terme in Sinus-Reihenentwlickung
N = 25;

% Anzahl Schleifendurchläufe
M = 2000;

% Ergebnisse
values = zeros(N * M, 1);

for l = 1:N
    for k = 1:M
        % gleichverteilte Sinuskoeffizienten generieren und auf [- 0.5, 0.5]
        % skalieren
        ucoeffs = 2 * (rand(l, 1) - 0.5);

        % H^{1}-Norm berechnen und normierte Funktion als Sinus-
        % Reihenentwicklung erzeugen
        unorm = h1_norm_sine(ucoeffs);
        u     = @(x) generate_sine_series(x, ucoeffs / unorm);

        % L_{\infty}-Norm bestimmen
        values((l - 1) * M + k) = linfty_norm(u);
    end
end

% Approximation an Einbettungskonstante ausgeben
sprintf('Einbettungskonstante: %f', max(values))
y = max(values);

% figure();
% scatter(1:length(values), values);
% title('Vergleich Abschätzung und experimentelle Stetigkeitskonstante');
% xlabel('l, wobei \phi_l(x) = c_l * sin(l \pi x)');
% legend('Abschätzung', 'Experimentell');
