function [B, ufun] = start( num_M, num_Q )
%% Galerkin-Verfahren f�r eindimensionale Raum-Zeit-Variationsformulierung.
% Das Variationsproblem $b(u, v) = f(v)$ mit $u \in \mathcal X$ und $v \in \mathcal Y$ wird durch
% Galerkin-Projektion auf endlichdimensionale Unterr�ume $\mathcal X_N$ und $\mathcal Y_N$
% approximiert.

%% Vorbereitung

% Varianten der Auswertung der einzelnen Integrale der Bilinearform
% Zur Auswahl stehen 0 und 1 f�r numerische Quadratur, respektive exakte Auswertung.
variante_term1 = 1;
variante_term3 = 1;

% Zeitintervall $I$
tspan = [0 1];
% Ortsintervall $\Omega$
xspan = [0 1];

% Anzahl der Sinus-Basisfunktionen f�r die Raumvariable $x$.
% num_M = 20;
% Anzahl der Legendre-Polynome f�r die Zeitvariable $t$.
% num_Q = 20;

% Gesamtzahl der Basisfunktionen f�r $\mathcal X_N$ festlegen.
num_X_j = num_M;
num_X_k = num_Q + 1;

% Gesamtzahl der Basisfunktionen f�r $\mathcal Y_N$ festlegen.
num_Y_l = num_M;
num_Y_m = num_Q;
num_Y_n = num_M;

% Die Dimensionen von $\mathcal X_N$ und $\mathcal Y_N$ m�ssen gleich sein um eine quadratische
% Steifigkeitsmatrix zu erhalten. (W�nschenswert!)
dim_X = num_X_j * num_X_k;
dim_Y = num_Y_l * num_Y_m + num_Y_n;
assert(dim_X == dim_Y)

%% Parameter der parabolischen PDE
% PDE lautet: u'(t) - c_0 \Delta u(t) + w u(t) = g(t), u(0) = u_0.
% Faktor vor Diffusionsterm
c_D = 1;

% Parameter f�r die Reihenentwicklung von \omega
epsilon = 1 / 100;
kappa = 99 / 100;

% Riemann-Zeta(2+epsilon), approx.
zeta = 1.65;

% Skalierung des Reaktionsterms
c_R = 1;

% Koeffizienten des Reaktionsterms
N_sigmas = 5;
% sigmas = c_R * (rand(N_sigmas + 1, 1) - 0.5);
sigmas = c_R * (ones(N_sigmas + 1, 1) - 0.5) ./ ((N_sigmas + 1):-1:1)';

% Check, ob hinreichende Bedingung f�r Regularit�t in sigma vorliegt.
assert(c_R < (c_D * kappa * pi^(4 + epsilon)) / (pi^(2 + epsilon) + 4));

% "Parameter-Funktion" $\omega$ festlegen.
omega = @(x) omega_sinus(x, 1/100, N_sigmas, sigmas);

% grid = linspace(0, 1, 100);
% plot(grid, w(grid));
% return;

% Anfangsbedingung $u0$ festlegen.
% u0 = @(x) 1 / 2 * sin(pi*x);
u0 = @(x) -(x) .* (x-1) .* (x+3);

% Und den Quellterm $g$ ebenfalls festlegen.
g = @(t, x) 0 + 0 .* t + 0 .* x;

%% Aufbau der Steifigkeitsmatrix und des Lastvektors

% Listen f�r die Erstellung der Sparse-Matrix $B$ erstellen.
% TODO: Eventuell noch Laufzeit rausholbar, indem diese mit einer (noch zu n�her zu bestimmenden
% Gr��e) preallokiert werden.
idx = [];
idy = [];
entry = [];

% �ber alle Basisfunktionen iterieren f�r $\mathcal X_N$ und $\mathcal Y_N$ iterieren. Dabei nutzen
% wir die Kreuzprodukt-Struktur von $\mathcal Y_N$ aus, indem wir die beiden Komponenten einzeln
% betrachten.
for j = 1:num_X_j
    for k = 1:num_X_k
        % Nur die erste Komponente von $v$, das hei�t $v = (v_1, 0)$.
        for l = 1:num_Y_l
            for m = 1:num_Y_m
                % Die Bilinearform $b(\cdot, \cdot)$ wurde in die einzelnen Terme zerlegt, so dass
                % m�glichst viel bereits per Hand vorberechnet werden kann.

                % $\int_T \int_\Omega u_t(t) v_1(t) dx dt$
                term1 = 0;
                if variante_term1 == 0
                    % Variante 1: Integral mittels numerischer Quadratur approximieren.
                    if j == l
                        term1 = integral(@(t) legendre_dP(t, k - 1) .* ...
                            legendre_P(t, m - 1), tspan(1), tspan(2)) / 2;
                    end;
                elseif variante_term1 == 1
                    % Variante 2: Ausnutzen, dass die Legendre-Polynome eine Orthogonalbasis f�r die
                    % Polynome darstellen. Damit vereinfacht sich alles zu...
                    if j == l
                        if k == 2 && m == 1
                            term1 = 1;
                        elseif k > m && mod(k, 2) ~= mod(m, 2)
                            term1 = 1 / 2;
                        end
                    end
                else
                    assert(false);
                end

                % $\int_T \int_\Omega c_D \nabla u(t) \nabla v_1(t) dx dt$
                term2 = 0;
                if j == l && k == m
                    term2 = c_D * ((pi * j)^2 / (2 * (k - 1) + 1)) / 2;
                end;

                % $\int_T \int_\Omega \omega u(t) v_1(t) dx dt$
                if variante_term3 == 0
                    % Variante 1: Integral mittels numerischer Quadratur approximieren.
                    term3 = 0;
                    if k == m
                        term3 = (1 / (2 * (k - 1) + 1)) * ...
                            integral(@(x) omega(x) .* sin(pi*j*x) .* ...
                            sin(pi*l*x), xspan(1), xspan(2));
                    end;
                elseif variante_term3 == 1
                    % Variante 2: Integral vorher explizit l�sen und hier nur noch auswerten.
                    % Problem: das Integral vom Produkt dreier Sinusfunktionen ist alles andere als
                    % sch�n. Dadurch ergeben sich hier einige Fallunterscheidungen...
                    term3 = 0;
                    if k == m
                        % Konstanter Term ist noch einfach
                        if j == l
                            term3 = sigmas(1) / 2;
                        end
                        % Die Sinusfunktionen dagegen nicht mehr. Hier wurde stark vereinfacht und
                        % alle F�lle, die automatisch Null ergeben direkt weggelassen
                        for iter = 1:N_sigmas
                            tmp = 0;
                            if mod(iter + j + l, 2) == 1
                                denom = (iter - j - l)*(iter + j - l)*(iter - j + l)*(iter + j + l) * pi;
                                if denom ~= 0
                                    tmp = (- 4 * iter * j * l) / denom;
                                end
                            end
                            term3 = term3 + (sigmas(iter + 1) / (iter * pi)^(1 + epsilon)) * tmp;
                        end
                    end;
                    term3 = (1 / (2 * (k - 1) + 1)) * term3;
                else
                    assert(false);
                end;

                % Der berechnete Wert wird nat�rlich nur dann abgespeichert, wenn er ungleich Null
                % ist.
                value = term1 + term2 + term3;
                if value ~= 0
                    x_pos = (j - 1) * num_X_k + k;
                    y_pos = (l - 1) * num_Y_m + m;

                    idx(end + 1) = x_pos;
                    idy(end + 1) = y_pos;
                    entry(end + 1) = value;
                end
            end
        end

        % Und jetzt die zweite Komponente von $v$, das hei�t $v = (0, v_2)$.
        for n = 1:num_Y_n
            if j == n
                tmp4 = (-1)^(k - 1) / 2;

                x_pos = (j - 1) * num_X_k + k;
                y_pos = num_Y_l * num_Y_m + n;

                idx(end + 1) = x_pos;
                idy(end + 1) = y_pos;
                entry(end + 1) = tmp4;
            end;
        end
    end
end
% Steifigkeitsmatrix aus den berechneten Werten erzeugen.
B = sparse(idy, idx, entry, dim_Y, dim_X);
return;

% F�r den Lastvektor erzeugen wir drei neue Hilfslisten
idx = [];
idy = [];
entry = [];

% Wie zuvor betrachten wir die beiden Komponenten von $\mathcal Y_N$ getrennt. Zun�chst also $v =
% (v_1, 0)$.
for l = 1:num_Y_l
    l
    for m = 1:num_Y_m
        value = integral2(@(t, x) g(t, x) .* sin(pi*l*x) .* ...
            legendre_P(t, m - 1), ...
            tspan(1), tspan(2), xspan(1), xspan(2));
        if value ~= 0
            idx(end + 1) = (l - 1) * num_Y_m + m;
            idy(end + 1) = 1;
            entry(end + 1) = value;
        end
    end
end
% Und nun $v = (0, v_2)$.
for n = 1:num_Y_n
    value = integral(@(x) u0(x) .* sin(pi*n*x), xspan(1), xspan(2));
    if value ~= 0
        idx(end + 1) = num_Y_l * num_Y_m + n;
        idy(end + 1) = 1;
        entry(end + 1) = value;
    end
end

F = sparse(idx, idy, entry, dim_Y, 1);

%% Das lineare Gleichungssystem $B u = F$ l�sen und die L�sung aufbauen.

% Dazu wird zun�chst das LGS gel�st.
u = B \ F;

% Und anschlie�end der daraus resultierende Koeffizientenvektor verwendet um die L�sung zu
% rekonstruieren.
ufun = @(t, x) reconstruct_solution(t, x, num_X_j, num_X_k, u);

% Und weil's so toll ist, plotten wir die L�sung bei Bedarf auch noch.
t_plot = 1;
if t_plot
    figure(2)
    tgrid = linspace(tspan(1), tspan(2), 100);
    xgrid = linspace(xspan(1), xspan(2), 50);
    [T, X] = meshgrid(tgrid, xgrid);
    mesh(T, X, ufun(T, X));
end
