function [MX, MY] = diskrete_normen(opt)
% DISKRETE_NORMEN Berechnet die Massematritzen für die diskreten Normen.

% X-Norm
% Erster Anteil: $L_2(0, T; H^1_0(\Omega))$
vecX1 = zeros(opt.dim.X, 1);
for j = 1:opt.num.Xj
    for k = 1:opt.num.Xk
        pos = (j - 1) * opt.num.Xk + k;
        vecX1(pos) = (1 + (pi * j)^2) / (2 * (2*(k - 1) + 1));
    end
end

% Zweiter Anteil: $L_2(0, T; H^{-1}(\Omega))$
blockX2 = zeros(opt.num.Xk);
for k1 = 1:opt.num.Xk
    for k2 = 1:opt.num.Xk
        if mod(k1 + k2, 2) == 0
            blockX2(k1, k2) = legendre_dP(1, k2 - 1);
            blockX2(k2, k1) = legendre_dP(1, k2 - 1);
        end
    end
end

% Aufaddieren und gut...
MX = spdiags(vecX1, 0, opt.dim.X, opt.dim.X) + ...
		kron(speye(opt.num.Xj), blockX2);

% Y-Norm
% Erster Anteil: $L_2(0, T; H^1_0(\Omega))$
vecY = zeros(opt.dim.Y, 1);
for l = 1:opt.num.Yl
    for m = 1:opt.num.Ym
        pos = (l - 1) * opt.num.Ym + m;
        vecY(pos) = (1 + (pi * l)^2) / (2 * (2*(m - 1) + 1));
    end
end

% Zweiter Anteil: $L_2(\Omega)$
vecY((opt.num.Yl * opt.num.Ym + 1):end) = 1 / 2;

% und gut!
MY = spdiags(vecY, 0, opt.dim.Y, opt.dim.Y);

end
