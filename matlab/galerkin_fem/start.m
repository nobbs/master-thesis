data = {};
data.xspan = [0 1];
data.tspan = [0 1];

dim = {};
dim.K = 1099;
dim.J = 1098;
dim.M = dim.J;
dim.N = dim.K * dim.J  + dim.M;

dt = (data.tspan(2) - data.tspan(1)) / dim.K;

% FIXME: remove
% [T1, T2] = time(dim);
% B = kron(N, V) + kron(M, A);

% Assembly
M = spacev(data.xspan, dim.M);
[A, A1, A2] = bfa(data.xspan, dim.J);

u = zeros(dim.M, dim.K + 1);

% u0 = @(x) sin(pi * x);
% u0 = @(x) -(x) .* (x-1) .* (x+3);
t = u0(linspace(data.xspan(1), data.xspan(2), dim.M+2));

% Startwert
u(:, 1) = t(2:end-1);

% rechte Seite
g = zeros(dim.M, 1);

M1 = (M / dt + A / 2);
M2 = (M / dt - A / 2);

% Crank-Nicolson...
for k = 2:(dim.K + 1),
	up = u(:, k - 1);
	uc = M1 \ (g + M2 * up);
	u(:, k) = uc;
end

size(u)

% Homogene Randbedingung dranpadden
upad = zeros(dim.M + 2, dim.K + 1);
upad(2:end-1, :) = u;

% Plotten
tgrid = linspace(data.tspan(1), data.tspan(2), dim.K + 1);
xgrid = linspace(data.xspan(1), data.xspan(2), dim.M + 2);
% 
[T, X] = meshgrid(tgrid, xgrid);
% mesh(T, X, upad);
