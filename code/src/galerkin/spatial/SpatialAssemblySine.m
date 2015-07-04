classdef SpatialAssemblySine < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a sine basis. The
  % resulting spatial discretization has homogenuous boundary conditions.
  %
  % @todo test me

  properties
    % witdh of the spatial interval @type double @type double
    xwidth;
  end

  methods

    function obj = SpatialAssemblySine(xwidth)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   xwidth: width of the spatial interval @type double

      if nargin > 0
        obj.xwidth = xwidth;
      end
    end

    function Mx = massMatrix(obj, nX, nY)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`
      % where `\sigma` are sine basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      Mx    = spdiags(obj.xwidth * ones(min(nX, nY), 1) / 2, 0, nY, nX);
    end

    function Ax = stiffnessMatrix(obj, nX, nY)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l` where `\sigma` are sine basis functions.
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
      for jdx = 1:min(nX, nY)
        tmp(jdx)  = (pi * jdx)^2 / (2 * obj.xspan(2));
      end
      Ax = spdiags(tmp, 0, nY, nX);
    end

  end % methods

end % classdef
