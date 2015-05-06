function [M, N] = time(dim)
	% Zeitspanne
	tspan = [0, 1];
	dt = (tspan(2) - tspan(1)) / dim.K;

	% sigma und tau
	M = zeros(dim.K, dim.K);
	% sigma' und tau
	N = zeros(dim.K, dim.K);

	for k = 1:dim.K
		M(k, k) = dt / 2;
		if k < dim.K
			M(k, k + 1) = dt / 2;
		end

		N(k, k) = 1;
		if k < dim.K
			N(k, k + 1) = -1;
		end
	end
end
