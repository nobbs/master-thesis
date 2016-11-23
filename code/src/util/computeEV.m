function [lm, lmvec, retries] = computeEV(A, B, shift)
  % Compute the largest generalized eigenvalue.
  %
  % This method computes the largest eigenvalue and the respective
  % eigenvector of the generalized eigenvalue problem `(\mat{A} +
  % \text{shift} * \mat{B})\vec{x} = \lambda \mat{B} \vec{x}`.
  %
  % Attention:
  %   There is some error handling included for the case that the matlab
  %   function eigs can't get the eigenvalue with default settings. So it
  %   may take some time, but it should get the eigenvalue most of the time,
  %   but possibly with a bad accuracy.
  %
  % Parameters:
  %   A: matrix on the left side of the eigenvalue problem @type matrix
  %   B: matrix on the right side of the eigenvalue problem @type matrix
  %   shift: shifts A on the left side by shift * B. this is needed if you
  %     want to compute the smallest eigenvalue; in this case you have to
  %     set shift to the value of the largest eigenvalue. @type double
  %
  % Return values:
  %   lm: largest eigenvalue @type double
  %   lmvec: eigenvector of the largest eigenvalue @type colvec
  %
  % @todo improve the error handling / accuracy decrease options

  retries = 0;

  % change the `eigenvalue won't converge` warning to an error, so that we
  % can catch it.
  origstate = warning;
  warning('error', 'MATLAB:eigs:NoEigsConverged'); %#ok<CTPCT>

  % setup the options structure
  opts = struct();
  opts.tol = eps;

  try
    % we first try to solve the eigenvalue problem with default settings.
    [lmvec, lm] = eigs(A - shift * B, B, 1, 'lm');
    % revert the shift
    lm = lm + shift;
  catch err
    % if the the default settings won't work, we have to try again with some
    % modifications
    if strcmp(err.identifier, 'MATLAB:eigs:NoEigsConverged') || strcmp(err.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
      % flag to check whether we found a eigenvalue
      hasEV = false;
      % decrease the accuracy
      opts.tol = opts.tol * 100;
      % set the number of used Lanczos vectors
      opts.p = 50;
      while ~hasEV
        try
          retries = retries + 1;
          % and now try again
          [lmvec, lm] = eigs(A - shift * B, B, 1, 'lm', opts);
          lm = lm  + shift;
          hasEV = true;
        catch err2
          hasEV = false;
          if strcmp(err2.identifier, 'MATLAB:eigs:NoEigsConverged') || strcmp(err2.identifier, 'MATLAB:eigs:ARPACKroutineErrorMinus14')
            % decrease the accuracy
            opts.tol = opts.tol * 100;
            % set the number of used Lanczos vectors
            opts.p = 50;
          else
            % other errors should be passed back
            rethrow(err2);
          end
        end
      end
    else
      % other errors should be passed back
      rethrow(err);
    end
  end

  % if we found an eigenvalue that looks like it's zero, it probably is
  % zero, so ...
  if abs(lm) < sqrt(eps)
    lm = 0;
  end

  % revert the warning state changes
  warning(origstate);
end
