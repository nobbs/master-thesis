function [lm, lmvec, retries] = computeMaxEV(A, B)
  [lmvec, lm] = eigs(A, B, 1, 'lm');
  retries = 0;
end