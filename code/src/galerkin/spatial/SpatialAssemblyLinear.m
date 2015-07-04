classdef SpatialAssemblyLinear < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a piecewise linear
  % nodal basis. The given methods support homogenuous and periodic boundary
  % conditions.
  %
  % @todo fixme

  methods(Static)

    function Mx = massMatrix(nK, xgrid, bc)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`,
      % which are piecewise linear nodal basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %   xgrid: grid of the spatial interval @type vector
      %   bc: type of boundary condition. possible values are 'homogenuous' and
      %     'periodic'. @type string @default 'homogenuous'
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix
      %
      % @todo fixme

      if nargin == 2
        bc = 'homogenuous';
      end

      switch bc
        case 'homogenuous'
          % @todo fixme
          D   = diff(xgrid);
          Mx = spdiags([D 0; 0 D]' * [1/6 1/3 0; 0 1/3 1/6], [-1 0 1], nK, nK);
          % modify the matrix to add the homogenuous boundary condition
          Mx(1, 2) = 0;
          Mx(1, 1) = 1;
          Mx(end, end - 1) = 0;
          Mx(end, end) = 1;
        case 'periodic'
          % @todo fixme
          tmp    = xwidth * ones(min(nX, nY), 1) / 2;
          tmp(1) = tmp(1) * 2;
          Mx    = spdiags(tmp, 0, nY, nX);
        otherwise
          error();
      end
    end

    function Ax = stiffnessMatrix(nK, xgrid, bc)
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
      % @todo fixme

      iD  = 1 ./ diff(xgrid);
      Ax = spdiags([iD 0; 0 iD]' * [-1 1 0; 0 1 -1], [-1 0 1], nK, nK);
    end

  end

end % classdef
