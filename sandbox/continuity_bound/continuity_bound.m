% Berechnet die Stetigkeits-Konstante für $A_j$ mit 
% $\phi_j = \frac{\sin(\pi * j * x)}{(j*\pi)^{1 + \epsilon})}.

% Anzahl Sinuskoeffizienten (insgesamt 2 * N + 1)
N = 25;

% Anzahl Schleifendurchläufe
M = 2500;

% Zu testendes $\phi_j$
l       = 20;
epsilon = 1 / 100;
K       = 5;
phi     = @(x) K * sin(l*pi*x) / (l*pi)^(1+epsilon);

xspan = [0, 1];
values = zeros(M, 1);

for k = 1:M
    if mod(k, 50) == 0
        k
    end
    
    % zufällig Sinuskoeffizienten generieren
    ucoeffs = rand(N, 1);
    vcoeffs = rand(N, 1);
    
    % Skalarprodukt, Normen und Schranke bestimmen
    unorm = h1_norm_sine(ucoeffs);
    vnorm = h1_norm_sine(vcoeffs);

    % Funktionen aus Sinuskoeffizienten erzeugen
    u  = @(x) generate_sine_series(x, ucoeffs / unorm);
    v  = @(x) generate_sine_series(x, vcoeffs / vnorm);
    
    values(k) = abs(l2_skp(@(x) phi(x) .* u(x), v, xspan));
end

figure()
max(values)
scatter(1:M, values);
