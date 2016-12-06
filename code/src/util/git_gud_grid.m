function points = git_gud_grid( pspan, n, dim )

N = floor(n^(1/dim));
x = linspace(pspan(1), pspan(2), N);

if dim == 2
  [X1, X2] = ndgrid(x);
  points = zeros(N^dim, dim);
  points(:, 1) = reshape(X1, [], 1);
  points(:, 2) = reshape(X2, [], 1);
elseif dim == 3
  [X1, X2, X3] = ndgrid(x);
  points = zeros(N^dim, dim);
  points(:, 1) = reshape(X1, [], 1);
  points(:, 2) = reshape(X2, [], 1);
  points(:, 3) = reshape(X3, [], 1);
elseif dim == 4
  [X1, X2, X3, X4] = ndgrid(x);
  points = zeros(N^dim, dim);
  points(:, 1) = reshape(X1, [], 1);
  points(:, 2) = reshape(X2, [], 1);
  points(:, 3) = reshape(X3, [], 1);
  points(:, 4) = reshape(X4, [], 1);
end

end

