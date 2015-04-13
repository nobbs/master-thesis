function [data, opt, omega] = setup(num_M, num_Q)
% SETUP Einstellungen und Bedingungen

%% Vorbereitung
opt = {};
data = {};

%% Optionen
% Varianten der Auswertung der einzelnen Integrale der Bilinearform
% Zur Auswahl stehen 0 und 1 für numerische Quadratur, respektive exakte Auswertung.
opt.variante_term1 = 1;
opt.variante_term3 = 1;

% Anzahl der Sinus-Basisfunktionen für die Raumvariable $x$.
opt.num.M = num_M;
% Anzahl der Legendre-Polynome für die Zeitvariable $t$.
opt.num.Q = num_Q;

% Gesamtzahl der Basisfunktionen für $\mathcal X_N$ festlegen.
opt.num.Xj = opt.num.M;
opt.num.Xk = opt.num.Q + 1;

% Gesamtzahl der Basisfunktionen für $\mathcal Y_N$ festlegen.
opt.num.Yl = opt.num.M;
opt.num.Ym = opt.num.Q;
opt.num.Yn = opt.num.M;

% Die Dimensionen von $\mathcal X_N$ und $\mathcal Y_N$ müssen gleich sein um eine quadratische
% Steifigkeitsmatrix zu erhalten. (Wünschenswert!)
opt.dim.X = opt.num.Xj * opt.num.Xk;
opt.dim.Y = opt.num.Yl * opt.num.Ym + opt.num.Yn;
assert(opt.dim.X == opt.dim.Y)

%% Daten
% Zeitintervall $I$
data.tspan = [0 1];
% Ortsintervall $\Omega$
data.xspan = [0 1];

%% Parameter der parabolischen PDE
% PDE lautet: u'(t) - c_0 \Delta u(t) + w u(t) = g(t), u(0) = u_0.
% Faktor vor Diffusionsterm
data.c_D = 1;
% Skalierung des Reaktionsterms
data.c_R = 1;

% Parameter für die Reihenentwicklung von \omega
data.eps = 1 / 100;
data.kappa = 99 / 100;

% Riemann-data.Zeta(2+data.eps), approx.
data.zeta = 1.65;

% Koeffizienten des Reaktionsterms
omega = {};
omega.N = 5;
% omega.sigmas = data.c_R * (rand(omega.N + 1, 1) - 0.5);
omega.sigmas = data.c_R * (ones(omega.N + 1, 1) - 0.5) ./ ((omega.N + 1):-1:1)';
% Check, ob hinreichende Bedingung für Regularität in sigma vorliegt.
assert(data.c_R < (data.c_D * data.kappa * pi^(4 + data.eps)) / (pi^(2 + data.eps) + 4));
% "Parameter-Funktion" $\omega$ festlegen.
omega.fun = @(x) omega_sinus(x, 1/100, omega.N, omega.sigmas);

% Anfangsbedingung $u0$ festlegen.
data.u0 = @(x) -(x) .* (x-1) .* (x+3);

% Und den Quellterm $g$ ebenfalls festlegen.
data.g = @(t, x) 0 + 0 .* t + 0 .* x;

end
