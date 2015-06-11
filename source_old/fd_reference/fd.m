function [Y, Tg, Xg] = fd(Nx, kind, omegafun)
%FD - Finite-Differenzen Implementierung für verschiedene Varianten der PPDE.
% Parameter:
%
% - Nx:
%   Anzahl Gitterpunkte für die Ortsdiskretisierung
% - kind:
%   Art der Randbedingung. Möglich: 'homogen', 'periodic'.
%   FIXME: Not Implemented!

tspan = [0 1];
xspan = [0 2];

% Steifigkeitsmatrix und Lastvektor erzeugen
[B1, B2, F] = assemble(Nx, kind, omegafun, xspan);

switch kind
  case 'homogen'
    % Im Falle homogener Randbedingungen werden die Randpunkte nicht als
    % Teil des Gleichungssystems bestimmt, sondern festgelegt. Dies
    % erfordert ein nachträgliches Padding der Lösung mit den Null-
    % Randwerten.

    Nxh = Nx - 2;
    B = B1 + B2;

    % Anfangswerte
    xgrid = linspace(xspan(1), xspan(2), Nx);
    xwidth = xgrid(2) - xgrid(1);
    u0 = sin(2 * pi * xgrid / 2);
    u0 = u0(2:end-1);
    % u0 = ones(Nxh, 1);

    % Verwende ode23s für die Lösung, alternativ wäre auch ein simples Timestepping im Style des expliziten Euler-Verfahrens denkbar.
    [T, Uh] = ode23s(@(t, u) B * u - F, tspan, u0);

    U = zeros(size(Uh, 1), size(Uh, 2) + 2);
    U(:, 2:end-1) = Uh;

    [Tg, Xg] = meshgrid(T, xgrid);

    % FIXME: Kontrollplot
    mesh(Tg, Xg, U');
    Y = U'
  case 'periodic'
    % Im Falle periodischer Randbedingungen ist die Lösung nur bis auf
    % Addition mit einer Konstante eindeutig und muss daher durch eine
    % zusätzliche Bedingung versehen werden um eine eindeutige Lösung zu
    % erhalten, beispielsweise das der Mittelwert gleich Null sei.
    Nxp = Nx - 1;
    B = B1 + B2;

    % Anfangswerte
    xgrid = linspace(xspan(1), xspan(2), Nx);
    xwidth = xgrid(2) - xgrid(1);
    % u0 = sin(2 * pi * xgrid(1:end-1));
    % u0 = u0(2:end-1);
    u0 = ones(Nxp, 1);

    % Verwende ode23s für die Lösung, alternativ wäre auch ein simples Timestepping im Style des expliziten Euler-Verfahrens denkbar.
    [T, Up] = ode23s(@(t, u) B * u - F, tspan, u0);

    U = zeros(size(Up, 1), size(Up, 2) + 1);
    U(:, 1:end-1) = Up;
    U(:, end) = Up(:, 1);
    % FIXME: Kontrollplot
    mesh(U);
  otherwise
    error('Not implemented!');
end

end

function [B1, B2, F] = assemble(Nx, kind, omegafun, xspan)
%ASSEMBLE - Erzeugt Steifigkeitsmatrix und Lastvektor mittels Finiter-
% Differenzen.
%
% Parameter:
% - Nx:
%   Anzahl Gitterpunkte für die Ortsdiskretisierung
% - kind:
%   Art der Randbedingung. Möglich: 'homogen', 'periodic'.
%   FIXME: Not Implemented!

switch kind
  case 'homogen'
    Nxh = Nx - 2;
    xgrid = linspace(xspan(1), xspan(2), Nx);
    xwidth = xgrid(2) - xgrid(1);

    % Steifigkeitsmatrix aufstellen
    % Erster Anteil: Laplace-Operator mit Null-Randbedingung
    B1 = spdiags(repmat([1, -2, 1], Nxh), [-1, 0, 1], Nxh, Nxh) / xwidth^2;

    % Zweiter Anteil: "Reaktion"
    % B2 = assemble_omega_from_vector(Nx, kind, omegafun(xgrid(2:end-1)));
    B2 = sparse(1,1);

    % Lastvektor aufstellen
    F = zeros(Nxh, 1);
  case 'periodic'
    Nxp = Nx - 1;
    xgrid = linspace(xspan(1), xspan(2), Nx);
    xwidth = xgrid(2) - xgrid(1);

    % Steifigkeitsmatrix aufstellen
    % Erster Anteil: Laplace-Operator mit periodischen Randbedingungen
    B1 = spdiags(repmat([1, -2, 1], Nx), [-1, 0, 1], Nxp, Nxp) / xwidth^2;
    B1(1, Nxp) = 1 / xwidth^2;
    B1(Nxp, 1) = 1 / xwidth^2;

    % Zweiter Anteil: "Reaktion"
    B2 = assemble_omega_from_vector(Nx, kind, omegafun(xgrid(1:end-1)));

    % Lastvektor aufstellen
    F = zeros(Nxp, 1);
  otherwise
    error('Not implemented!');
end

end

function Bomega = assemble_omega_from_vector(Nx, kind, omega)
switch kind
  case 'homogen'
    Nxh = Nx - 2;
    Bomega = - spdiags(omega, 0, Nxh, Nxh);
  case 'periodic'
    Nxp = Nx - 1;
    Bomega = - spdiags(omega, 0, Nxp, Nxp);
  otherwise
    error('Not implemented!');
end
end

function Bomega = assemble_omega_from_series_coeffs(Nx, kind, omega, kind_of_series)
error('Not Implemented!');
end
