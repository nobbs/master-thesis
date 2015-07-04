classdef SpatialAssemblyFourier < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a fourier basis. The
  % resulting spatial discretization has periodic boundary conditions.

  properties
    % witdh of the spatial interval @type double @type double
    xwidth;
  end

  methods

    function obj = SpatialAssemblyFourier(xwidth)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   xwidth: width of the spatial interval @type double

      obj.xwidth = xwidth;
    end

    function Mx = massMatrix(obj, nX, nY)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`
      % where `\sigma` are Fourier basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      tmp    = obj.xwidth * ones(min(nX, nY), 1) / 2;
      tmp(1) = tmp(1) * 2;
      Mx     = spdiags(tmp, 0, nY, nX);
    end

    function Ax = stiffnessMatrix(obj, nX, nY)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l` where `\sigma` are Fourier basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   Ax: spatial stiffness matrix @type matrix
      %
      % @todo optimize

      tmp = zeros(min(nX, nY), 1);
      for jdx = 2:min(nX, nY)
        if mod(jdx, 2) == 0
          tmp(jdx) = (pi * jdx)^2 / ( 2 * obj.xwidth);
        else
          tmp(jdx) = (pi * (jdx - 1))^2 / ( 2 * obj.xwidth);
        end
      end
      Ax = spdiags(tmp, 0, nY, nX);
    end

  end % methods

end % classdef
