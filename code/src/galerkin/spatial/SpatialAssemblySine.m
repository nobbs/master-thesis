classdef SpatialAssemblySine < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a sine basis. The
  % resulting spatial discretization has homogenuous boundary conditions.
  %
  % @todo test me

  methods(Static)

    function Mx = massMatrix(nX, nY, xwidth)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`,
      % which are sine basis functions in this case.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   nY: the limit of iteration for the `l` index @type integer
      %   xwidth: width of the spatial interval @type double
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix

      Mx    = spdiags(xwidth * ones(min(nX, nY), 1) / 2, 0, nY, nX);
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
      for jdx = 1:min(nX, nY)
        tmp(jdx)  = (pi * jdx)^2 / (2 * obj.xspan(2));
      end
      Ax = spdiags(tmp, 0, nY, nX);
    end

  end % methods

end % classdef
