function [lm, lmvec, retries] = computeEVdense(A, B, shift)
  [V, D] = eig(full(A - shift * B), full(B));
  D = diag(D);
  [lm, i] = max(abs(D));
  lmvec = V(:, i);

  % if we found an eigenvalue that looks like it's zero, it probably is
  % zero, so ...
  if abs(lm) < sqrt(eps)
    lm = 0;
  end

  lm = lm + shift;
end
