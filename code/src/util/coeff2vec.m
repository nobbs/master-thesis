%% coeff2vec: function description
function vector = coeff2vec(coeffs, xgrid, N)
	vector = zeros(size(xgrid, 1), size(xgrid, 2));
	for idx = 1:N
		if mod(idx, 2) == 1
		% odd: cosine
			vector = vector + coeffs(idx) * cos(pi * (idx + 1) * xgrid / 10);
		else
		% even: sine
			vector = vector + coeffs(idx) * sin(pi * idx * xgrid / 10);
		end
	end
end
