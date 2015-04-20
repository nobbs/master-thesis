data = {};
data.xspan = [0 1];
data.tspan = [0 1];

dim = {};
dim.K = 500;
dim.J = 500;
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
u0 = @(x) -(x) .* (x-1) .* (x+3);
t = u0(linspace(data.xspan(1), data.xspan(2), dim.M+2));

% Startwert
u(:, 1) = t(2:end-1);

% rechte Seite
g = zeros(dim.M, 1);

for k = 2:(dim.K + 1),
	up = u(:, k - 1);
	uc = (M / dt + A / 2) \ (g + (M / dt - A / 2) * up);
	u(:, k) = uc;
end

size(u)

% Homogene Randbedingung dranpadden
upad = zeros(dim.M + 2, dim.K + 1);
upad(2:end-1, :) = u;

% Plotten
T = linspace(data.tspan(1), data.tspan(2), dim.K + 1);
X = linspace(data.xspan(1), data.xspan(2), dim.M + 2);

[T, X] = meshgrid(T, X);
mesh(T, X, upad);
