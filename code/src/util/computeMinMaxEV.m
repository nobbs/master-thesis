function [mi, ma, retries] = computeMinMaxEV(A, B)
  % Computes the minimal and maximal generalized eigenvalues.
  %
  % Parameters:
  %   A: matrix on the left side @type matrix
  %   B: matrix on the right side @type matrix
  %
  % Return values:
  %   mi: minimal eigenvalue @type double
  %   ma: maximal eigenvalue @type double

  % try new stuff:
  [lm, ~, retries1] = computeMaxEV(A, B);
  if lm >= 0
    ma = lm;
    [mi, ~, retries2] = computeMinEV(A, B);
  else
    mi = lm;
    [ma, ~, retries2] = computeMinEV(A, B);
  end

  % first we calculate the eigenvalue with the largest absolute value
  % [lm, ~, retries1] = computeEV(A, B, 0);
  % % if it's positive, then it's the maximal eigenvalue, else the
  % % minimal eigenvalue
  % if lm >= 0
  %   ma = lm;
  %   % now shift to get the minimal eigenvalue
  %   [mi, ~, retries2] = computeEV(A, B, lm);
  % else
  %   mi = lm;
  %   % now shift to get the maximal eigenvalue
  %   [ma, ~, retries2] = computeEV(A, B, lm);
  % end

  retries = retries1 + retries2;
end