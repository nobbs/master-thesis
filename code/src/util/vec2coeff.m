%% vec2coeff: function description
function coeffs = vec2coeff(vector, xgrid, N)
	coeffs = zeros(N, 1);
	for idx = 1:N
		if mod(idx, 2) == 1
		% odd: cosine
			coeffs(idx) = trapez(vector .* cos(pi * (idx + 1) * xgrid / 10).', xgrid, 1) / 5;
		else
		% even: sine
			coeffs(idx) = trapez(vector .* sin(pi * idx * xgrid / 10).', xgrid, 1) / 5;
		end
	end
end
