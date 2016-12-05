function [lm, lmvec, retries] = computeMinEV(A, B)
  s = warning('error', 'MATLAB:nearlySingularMatrix');
  warning('error', 'MATLAB:eigs:SigmaNearExactEig');
  
  try
    [lmvec, lm] = eigs(A, B, 1, 'sm');
    retries = 0;
  catch
    lmvec = 0;
    lm = 0;
    retries = 0;
  end

  warning(s);
end