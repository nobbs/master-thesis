function [M] = spacev(xspan, M)
	dx = (xspan(2) - xspan(1)) / M

	% phi und phi
	M = spdiags(dx / 6 * repmat([1, 4, 1], M, 1), [-1 0 1], M, M);
end
