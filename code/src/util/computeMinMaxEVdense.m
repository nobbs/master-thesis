function [mi, ma, retries] = computeMinMaxEVdense(A, B)
  % Computes the minimal and maximal generalized eigenvalues.
  %
  % Parameters:
  %   A: matrix on the left side @type matrix
  %   B: matrix on the right side @type matrix
  %
  % Return values:
  %   mi: minimal eigenvalue @type double
  %   ma: maximal eigenvalue @type double

  D = eig(full(A), full(B));
  ma = max(abs(D));
  mi = min(abs(D));

  retries = 0;
end