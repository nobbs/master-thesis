classdef SpatialAssemblyFourier < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a fourier basis. The
  % resulting spatial discretization has periodic boundary conditions.

  methods(Static)

    function Mx = massMatrix(nX, nY, xwidth)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`,
      % which are Fourier basis functions in this case.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %   xwidth: width of the spatial interval @type double
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      tmp    = xwidth * ones(min(nX, nY), 1) / 2;
      tmp(1) = tmp(1) * 2;
      Mx     = spdiags(tmp, 0, nY, nX);
    end

    function Ax = stiffnessMatrix(nX, nY, xwidth)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l`.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %   xwidth: width of the spatial interval @type double
      %
      % Return values:
      %   Ax: spatial stiffness matrix @type matrix
      %
      % @todo optimize

      tmp = zeros(min(nX, nY), 1);
      for jdx = 2:min(nX, nY)
        if mod(jdx, 2) == 0
          tmp(jdx) = (pi * jdx)^2 / ( 2 * xwidth);
        else
          tmp(jdx) = (pi * (jdx - 1))^2 / ( 2 * xwidth);
        end
      end
      Ax = spdiags(tmp, 0, nY, nX);
    end

  end % methods

end % classdef
