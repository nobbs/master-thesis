classdef AssemblyGlobalAbstract < AssemblyAbstract
  % Common superclass for the assembly classes of global galerkin methods.
  %
  % The main purpose  of this class is to provide common shared methods which
  % are the same for the different types of galerkin methods with global spatial
  % and temporal basis function.

  properties
    % Number of spatial basis functions for the ansatz subspace @type integer
    nAnsatzSpatial;

    % Number of temporal basis polynomials for the ansatz subspace @type integer
    nAnsatzTemporal;

    % Number of spatial basis functions for the test subspace @type integer
    nTestSpatial;

    % Number of temporal basis polynomials for the test subspace @type integer
    nTestTemporal;

    % Number of spatial basis functions for the initial condition part of the
    % test subspace @type integer
    nTestSpatialIC;
  end

  properties(Dependent)
    % Dimension of the ansatz subspace @type integer
    dAnsatz;

    % Dimension of the test subspace @type integer
    dTest;
  end

  methods(Abstract)
    % Spatial basis functions.
    %
    % Evaluates the spatial basis function with the given index for the given
    % values of x. Can be used to define function handles and for numerical
    % integration.
    %
    % Parameters:
    %   index: index of the basis function @type integer
    %   x: values in which the function should be evaluated @type matrix
    %
    % Return values:
    %   val: values of the basis function in x @type matrix
    val = spatialBasisFunc(obj, index, x);

    % Temporal basis functions.
    %
    % Evaluates the temporal basis function with the given index for the given
    % values of t. Can be used to define function handles and for numerical
    % integration.
    %
    % Parameters:
    %   index: index of the basis function @type integer
    %   t: values in which the function should be evaluated @type matrix
    %
    % Return values:
    %   val: values of the basis function in t @type matrix
    val = temporalBasisFunc(obj, index, t);
  end

  methods

    function dim = get.dAnsatz(obj)
      dim = obj.nAnsatzSpatial * obj.nAnsatzTemporal;
    end

    function dim = get.dTest(obj)
      dim = obj.nTestSpatial * obj.nTestTemporal + obj.nTestSpatialIC;
    end

    function setNumberOfAnsatzFuncs(obj, nAnsatzSpatial, nAnsatzTemporal)
      % Set the number of basis functions for the ansatz subspace.
      %
      % Parameters:
      %   nAnsatzSpatial: Number of sine basis functions @type integer
      %   nAnsatzTemporal: Number of Legendre basis polynomials @type integer

      obj.nAnsatzSpatial  = nAnsatzSpatial;
      obj.nAnsatzTemporal = nAnsatzTemporal;
    end

    function setNumberOfTestFuncs(obj, nTestSpatial, nTestTemporal, nTestSpatialIC)
      % Set the number of basis functions for the test subspace.
      %
      % Parameters:
      %   nTestSpatial: Number of sine basis functions  @type integer
      %   nTestTemporal: Number of Legendre basis polynomials @type integer
      %   nTestSpatialIC: Number of sine basis functions for the initial condition
      %     @type integer

      obj.nTestSpatial   = nTestSpatial;
      obj.nTestTemporal  = nTestTemporal;
      obj.nTestSpatialIC = nTestSpatialIC;
    end

    function setNumberOfTestFuncsFromAnsatzFuncs(obj)
      % Set the number of basis functions for the test subspace according to the
      % number of basis functions for the subspace.
      %
      % This guaranties a quadratic system of linear equations.

      obj.nTestSpatial   = obj.nAnsatzSpatial;
      obj.nTestTemporal  = obj.nAnsatzTemporal - 1;
      obj.nTestSpatialIC = obj.nAnsatzSpatial;
    end

    function val = solutionFuncFromCoeffs(obj, solutionCoeffs, t, x)
      % Construct a function handle of the solution.
      %
      % Use it to define a solution function in (t, x) or evaluate the solution
      % directly for given t and x grids.
      %
      % Parameters:
      %   solutionCoeffs: solution vector of the linear system @type colvec
      %   t: temporal variable @type vecotr
      %   x: spatial variable @type vector
      %
      % Return values:
      %   solfun: function handle of the solution function @type function_handle

      val = zeros(size(t, 1), size(t, 2));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = cell(obj.nAnsatzSpatial, 1);
      for jdx = 1:obj.nAnsatzSpatial
        spatialValues{jdx} = obj.spatialBasisFunc(jdx, x);
      end
      temporalValues = cell(obj.nAnsatzTemporal, 1);
      for kdx = 1:obj.nAnsatzTemporal
        temporalValues{kdx} = obj.temporalBasisFunc(kdx, t);
      end

      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
        % Get the right coefficient
        pos = (jdx - 1) * obj.nAnsatzTemporal + kdx;

        %   % @todo normalization constant
        %   % cnrm = sqrt((1 + (pi * j)^2) / (2 * (2*(k - 1) + 1)) + ...
        %             % legendre_dP(1, k - 1));

        % evaluate the corresponding basis functions
        val = val + solutionCoeffs(pos) * spatialValues{jdx} .* temporalValues{kdx};
        end
      end
    end

  end

end
