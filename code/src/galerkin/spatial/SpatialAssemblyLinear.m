classdef SpatialAssemblyLinear < SpatialAssemblyAbstract
  % Assemble the spatial mass and stiffness matrices for a piecewise linear
  % nodal basis. The given methods support homogenuous and periodic boundary
  % conditions.
  %
  % @todo fixme

  properties
    % type of boundary condition. @type string
    % possible values are:
    %   homogenuous
    %   periodic
    bc;

    % grid of the spatial interval @type vector
    xgrid;
  end

  methods

    function obj = SpatialAssemblyLinear(xgrid, bc)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   xgrid: width of the spatial interval @type double @type vector
      %   bc: @type of boundary condition. possible: 'homogenuous', 'periodic'
      %     @type string

      obj.xgrid = xgrid;
      obj.bc    = bc;
    end

    function Mx = massMatrix(obj, nK)
      % Assemble the spatial mass matrix, that means we evaluate the integral
      % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`,
      % which are piecewise linear nodal basis functions.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %
      % Return values:
      %   Mx: spatial mass matrix @type matrix
      %
      % @todo fixme

      switch obj.bc
        case 'homogenuous'
          % @todo fixme
          D   = diff(obj.xgrid);
          Mx = spdiags([D 0; 0 D]' * [1/6 1/3 0; 0 1/3 1/6], [-1 0 1], nK, nK);
          % modify the matrix to add the homogenuous boundary condition
          Mx(1, 2) = 0;
          Mx(1, 1) = 1;
          Mx(end, end - 1) = 0;
          Mx(end, end) = 1;
        case 'periodic'
          % @todo fixme
        otherwise
          error();
      end
    end

    function Ax = stiffnessMatrix(obj, nK)
      % Assemble the spatial stiffness matrix, that means we evaluate the
      % integral `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j`
      % and `l`.
      %
      % Parameters:
      %   nX: the limit of iteration for the `j` index @type integer
      %
      % Return values:
      %   Ax: spatial stiffness matrix @type matrix
      %
      % @todo fixme

      iD  = 1 ./ diff(obj.xgrid);
      Ax = spdiags([iD 0; 0 iD]' * [-1 1 0; 0 1 -1], [-1 0 1], nK, nK);
    end

  end

end % classdef
