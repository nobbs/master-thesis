[B, ~] = steifigkeitsmatrix_homogen(data, omega, opt, data.tspan);
[F] = lastvektor_homogen(data, opt);

U = B \ F;

L = data.xspan(2);

tgrid = linspace(data.tspan(1), data.tspan(2), 100);
xgrid = linspace(data.xspan(1), data.xspan(2), 100);
[T, X] = meshgrid(tgrid, xgrid);
T = Tfd;
X = Xfd;

% u = reassemble_solution(U);
u = zeros(size(T, 1), size(T, 2));

for jdx = 1:opt.num.Xj
    for kdx = 1:opt.num.Xk
    	% Zur aktuellen Ansatzfunktion zugehörigen Koeffizienten bestimmen
    	cidx = (jdx - 1) * opt.num.Xk + kdx;
    	kk = kdx - 1;
    	% Normierungsfaktor bestimmen
		% cnrm = sqrt((1 + (pi * jdx)^2) / (2 * (2*(kk) + 1)) + ...
		            % legendre_dP_shifted(1, kk, data.tspan));
		% Ansatzfunktion auswerten
        u = u + U(cidx) * sin(pi * jdx * X / L) .* legendre_P_shifted(T, kk, data.tspan);
    end
end

mesh(T, X, u);
xlabel('T')
ylabel('X')
zlabel('U')
