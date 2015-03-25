% Berechnet die Stetigkeits-Konstante für $A_j$ mit
% $\phi_j = \frac{\sin(\pi * j * x)}{(j*\pi)^{1 + \epsilon})}.

% Anzahl Sinuskoeffizienten (insgesamt 2 * N + 1)
N = 25;

% Anzahl Schleifendurchläufe
M = 2000;

xspan = [0, 1];
values = zeros(N * M, 1);

for l = 1:N
l
    for k = 1:M
        % zufällig Sinuskoeffizienten generieren
        ucoeffs = 2 * (rand(l, 1) - 0.5);
        
        % Skalarprodukt, Normen und Schranke bestimmen
        unorm = h1_norm_sine(ucoeffs);
        
        % Funktionen aus Sinuskoeffizienten erzeugen
        u  = @(x) generate_sine_series(x, ucoeffs / unorm);

        values((l - 1) * M + k) = linfty_norm(u, 250);
    end
end

max(values)

figure();
scatter(1:length(values), values);
% title('Vergleich Abschätzung und experimentelle Stetigkeitskonstante');
% xlabel('l, wobei \phi_l(x) = c_l * sin(l \pi x)');
% legend('Abschätzung', 'Experimentell');
