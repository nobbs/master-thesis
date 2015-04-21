function y = exfun(u0fun, T, X, K)
	if nargin == 2
		K = 25;
	end

	y = zeros(size(T, 1), size(T, 2));
	for k = 1:K
		c = 2 * integral(@(x) u0fun(x) .* sin(pi * k * x), 0, 1);
		y = y + c * sin(pi * k * X) .* exp(-((k * pi)^2 + 1) * T);
	end
end
