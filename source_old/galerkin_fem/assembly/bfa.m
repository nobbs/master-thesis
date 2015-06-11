function [A, A1, A2] = bfa(xspan, J)

	dx = (xspan(2) - xspan(1)) / J

	% Konsante
	c1 = 1;
	c2 = 1;

	% phi' und phi'
	A1 = c1 * spdiags((1 / dx) * repmat([-1, 2, -1], J, 1), [-1 0 1], J, J);

	% phi und phi
	A2 = c2 * spdiags(dx / 6 * repmat([1, 4, 1], J, 1), [-1 0 1], J, J);

	A = A1 + A2;
end
