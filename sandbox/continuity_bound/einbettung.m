function [ y ] = einbettung()
% Berechnet die Einbettungskonstante für die Einbettung $H^{1}_{0}(\Omega)
% \subset L_{\infty}(\Omega)$ anhand von Sinusfunktionen mit zufällig
% gleichverteilten Koeffizienten.

% Intervallgröße
L = 5;

% Anzahl Terme in Sinus-Reihenentwlickung
N = 5;

% Anzahl Schleifendurchläufe
M = 200;

% Ergebnisse
values = zeros(N * M, 1);

pl = [];

for l = 1:N
    for k = 1:M
        % gleichverteilte Sinuskoeffizienten generieren und auf [- 0.5, 0.5]
        % skalieren
        ucoeffs = 2 * (rand(l, 1) - 0.5);

        % H^{1}-Norm berechnen und normierte Funktion als Sinus-
        % Reihenentwicklung erzeugen
        unorm = h1_norm_sine(ucoeffs, L);
        u     = @(x) generate_sine_series(x, ucoeffs / unorm, L);
        
        K = length(ucoeffs);
        ducoeffs = pi * (1:K)' .* ucoeffs / L;
        
        pl(end+1) = h1_norm(u, @(x) generate_cosine_series(x, ducoeffs / unorm, L), [0 L]);

%         plot(u(linspace(0, L, 100)));
%         linfty_norm(u, 100, L)
%         drawnow
        
        % L_{\infty}-Norm bestimmen
        values((l - 1) * M + k) = linfty_norm(u, 100, L);
    end
end

% plot(pl)
% return

% Approximation an Einbettungskonstante ausgeben
sprintf('Einbettungskonstante: %f', max(values))
y = max(values);

figure();
scatter(1:length(values), values);
title('Vergleich Abschätzung und experimentelle Stetigkeitskonstante');
xlabel('l, wobei \phi_l(x) = c_l * sin(l \pi x)');
% legend('Abschätzung', 'Experimentell');
