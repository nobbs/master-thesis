% Berechnet die Stetigkeits-Konstante für $A_j$ mit
% $\phi_j = \frac{\sin(\pi * j * x)}{(j*\pi)^{1 + \epsilon})}.

% 
Lspan = 1:25;

% Anzahl Sinuskoeffizienten (insgesamt 2 * N + 1)
N = 25;

% Anzahl Schleifendurchläufe
M = 2000;

% Zu testendes $\phi_j$
epsilon = 1 / 100;
K       = 5;

xspan = [0, 1];

bounds = zeros(length(Lspan), 1);
approxs = zeros(length(Lspan), 1);

for idx_j = 1:length(Lspan)
    j      = Lspan(idx_j);
    phi    = @(x) (K / (j * pi)^(1+epsilon)) * sin(j * pi * x);
    values = zeros(N * M, 1);
    max_val = 0;

    for l = 1:N
        for k = 1:M
            % zufällig Sinuskoeffizienten generieren
            ucoeffs = 2 * (rand(l, 1) - 0.5);
            vcoeffs = 2 * (rand(l, 1) - 0.5);
            
            % Skalarprodukt, Normen und Schranke bestimmen
            unorm = h1_norm_sine(ucoeffs);
            vnorm = h1_norm_sine(vcoeffs);
            
            % Funktionen aus Sinuskoeffizienten erzeugen
            u  = @(x) generate_sine_series(x, ucoeffs / unorm);
            v  = @(x) generate_sine_series(x, vcoeffs / vnorm);
            
%             values((l - 1) * M + k) = abs(l2_skp(@(x) phi(x) .* u(x), v, xspan));
            val     = abs(l2_skp(@(x) phi(x) .* u(x), v, xspan));
            max_val = max(val, max_val);
        end
    end
   
    % Abschätzung
    bounds(idx_j)  = (4 * K) / (j * pi)^(2 + epsilon);
    approxs(idx_j) = max_val;
end

figure();
semilogy(Lspan, bounds, Lspan, approxs);
title('Vergleich Abschätzung und experimentelle Stetigkeitskonstante');
xlabel('l, wobei \phi_l(x) = c_l * sin(l \pi x)');
legend('Abschätzung', 'Experimentell');
