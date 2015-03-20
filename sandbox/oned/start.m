%% Galerkin-Verfahren f�r eindimensionale Raum-Zeit-Variationsformulierung.
% Das Variationsproblem $b(u, v) = f(v)$ mit $u \in \mathcal X$ und
% $v \in \mathcal Y$ wird durch Galerkin-Projektion auf endlichdimensionale
% Unterr�ume $\mathcal X_N$ und $\mathcal Y_N$ approximiert.

%% Vorbereitung

% Zeitintervall $I$
tspan = [0 1];
% Ortsintervall $\Omega$
xspan = [0 1];

% Anzahl der Sinus-Basisfunktionen f�r die Raumvariable $x$.
num_M = 20;
% Anzahl der Legendre-Polynome f�r die Zeitvariable $t$.
num_Q = 20;

% Gesamtzahl der Basisfunktionen f�r $\mathcal X_N$ festlegen.
num_X_j = num_M;
num_X_k = num_Q + 1;

% Gesamtzahl der Basisfunktionen f�r $\mathcal Y_N$ festlegen.
num_Y_l = num_M;
num_Y_m = num_Q;
num_Y_n = num_M;

% Die Dimensionen von $\mathcal X_N$ und $\mathcal Y_N$ m�ssen gleich sein
% um eine quadratische Steifigkeitsmatrix zu erhalten. (W�nschenswert!)
dim_X = num_X_j * num_X_k;
dim_Y = num_Y_l * num_Y_m + num_Y_n;
assert(dim_X == dim_Y)

% "Parameter-Funktion" $\omega$ festlegen.
w = @(x) 1 + 0 .* x;

% Anfangsbedingung $u0$ festlegen.
u0 = @(x) x .* sin(pi*x);

% Und den Quellterm $g$ ebenfalls festlegen.
g = @(t, x) 0 .* t + 0 .* x;

%% Aufbau der Steifigkeitsmatrix und des Lastvektors

% Listen f�r die Erstellung der Sparse-Matrix $B$ erstellen.
% TODO: Eventuell noch Laufzeit rausholbar, indem diese mit einer (noch zu
% n�her zu bestimmenden Gr��e) preallokiert werden.
idx = [];
idy = [];
entry = [];

% �ber alle Basisfunktionen iterieren f�r $\mathcal X_N$ und $\mathcal Y_N$
% iterieren. Dabei nutzen wir die Kreuzprodukt-Struktur von $\mathcal Y_N$
% aus, indem wir die beiden Komponenten einzeln betrachten.
for j = 1:num_X_j
    for k = 1:num_X_k
        % Nur die erste Komponente von $v$, das hei�t $v = (v_1, 0)$.
        for l = 1:num_Y_l
            for m = 1:num_Y_m
                % Die Bilinearform $b(\cdot, \cdot)$ wurde in die einzelnen
                % Terme zerlegt, so dass m�glichst viel bereits per Hand
                % vorberechnet werden kann.

                % $\int_T \int_\Omega u_t(t) v_1(t) dx dt$
                tmp1 = 0;
                if j == l
                    tmp1 = integral(@(t) legendre_dP(t, k - 1) .* ...
                            legendre_P(t, m - 1), tspan(1), tspan(2)) / 2;
                end;

                % $\int_T \int_\Omega \nabla u(t) \nabla v_1(t) dx dt$
                tmp2 = 0;
                if j == l && k == m
                    tmp2 = ((pi * j)^2 / (2 * (k - 1) + 1)) / 2;
                end;

                % $\int_T \int_\Omega \omega u(t) v_1(t) dx dt$
                tmp3 = 0;
                if k == m
                    tmp3 = (1 / (2 * (k - 1) + 1)) * ...
                            integral(@(x) w(x) .* sin(pi*j*x) .* ...
                                    sin(pi*l*x), xspan(1), xspan(2));
                end;

                % Der berechnete Wert wird nat�rlich nur dann
                % abgespeichert, wenn er ungleich Null ist.
                value = tmp1 + tmp2 + tmp3;
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

% F�r den Lastvektor erzeugen wir drei neue Hilfslisten
idx = [];
idy = [];
entry = [];

% Wie zuvor betrachten wir die beiden Komponenten von $\mathcal Y_N$
% getrennt. Zun�chst also $v = (v_1, 0)$.
for l = 1:num_Y_l
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

% Und anschlie�end der daraus resultierende Koeffizientenvektor verwendet
% um die L�sung zu rekonstruieren.
ufun = @(t, x) reconstruct_solution(t, x, num_X_j, num_X_k, u);

% Und weil's so toll ist, plotten wir die L�sung bei Bedarf auch noch.
t_plot = 1;
if t_plot
    figure(2)
    tgrid = linspace(tspan(1), tspan(2), 33);
    xgrid = linspace(xspan(1), xspan(2), 33);
    [T, X] = meshgrid(tgrid, xgrid);
    mesh(T, X, ufun(T, X));
end
