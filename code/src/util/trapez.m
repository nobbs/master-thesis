function Q = trapez(f, grid, is_periodic)
  % Cumulative trapezoidal quadrature.
  %
  % Evaluates the integral `\int_{a}^{b} f(x) \diff x` through a cummulative
  % trapezoidal quadrature on a equidistant grid. If `f` is a periodic function
  % and the interval of integration `[a, b]` is a multiple of the period of `f`,
  % then is_periodic should be set to true and grid and function values should
  % be given for the half open interval `[a, b)`.
  %
  % Parameters:
  %   f: function values at the grid points @type vector
  %   grid: grid of the interval `[a, b]` @type vector
  %   is_periodic: toggles whether f is periodic on the given grid
  %     @default true @type logical
  %
  % Return values:
  %   Q: approximation of the integral value

  % set default values
  if nargin == 2
    is_periodic = 1;
  end

  % check whether vector of function values and grid have the same size
  assert(length(f) == length(grid));

  % compute the cumulative trapezoidal quadrature
  h = grid(2) - grid(1);
  if is_periodic
    Q = h * (f(1) + sum(f(2:end)));
  else
    Q = (h / 2) * (f(1) + f(end) + 2 * sum(f(2:end-1)));
  end

end
