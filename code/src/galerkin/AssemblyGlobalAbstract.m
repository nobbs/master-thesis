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

    % Normalize ansatz and test basis functions in the respective norms @type
    % logical
    useNormalization = true;

    % Store the norms of the non-normalized ansatz basis functions in the
    % respective norm @type vector
    AnsatzNormDiag;

    % Store the norms of the non-normalized test basis functions in the
    % respective norm @type vector
    TestNormDiag;
  end

  properties(Access = 'public')
    % holds the normalization of the ansatz functions in matrix form
    % @type sparsematrix
    ansatzNormalizationMatrix;

    % holds the normalization of the test functions in matrix form
    % @type sparsematrix
    testNormalizationMatrix;
  end

  properties(Dependent)
    % Dimension of the ansatz subspace @type integer
    nAnsatzDim;

    % Dimension of the test subspace @type integer
    nTestDim;
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

    % First derivative of spatial basis functions.
    %
    % Evaluates the first derivative of the spatial basis function with the
    % given index for the given values of x. Can be used to define function
    % handles and for numerical integration.
    %
    % Parameters:
    %   index: index of the basis function @type integer
    %   x: values in which the function should be evaluated @type matrix
    %
    % Return values:
    %   val: values of the basis function in x @type matrix
    val = spatialBasisFuncDerivative(obj, index, x)

    % Temporal basis functions.
    %
    % Evaluates the temporal basis function with the given index for the given
    % values of t. Can be used to define function handles and for numerical
    % integration.
    %
    % Parameters:
    %   index: index of the basis function @type integer
    %   t: values in which the function should be evaluated @type matrix
    %   tspan: custom temporal interval @type vector @default obj.tspan
    %
    % Return values:
    %   val: values of the basis function in t @type matrix
    val = temporalBasisFunc(obj, index, t, tspan);

    % First derivative of temporal basis functions.
    %
    % Evaluates the first derivative of a temporal basis function with the given
    % index for the given values of t. Can be used to define function handles
    % and for numerical integration.
    %
    % Parameters:
    %   index: index of the basis function @type integer
    %   t: values in which the function should be evaluated @type matrix
    %   tspan: custom temporal interval @type vector @default obj.tspan
    %
    % Return values:
    %   val: values of the basis function in t @type matrix
    val = temporalBasisFuncDerivative(obj, index, t, tspan)

    % Precompute the normalization constants of the ansatz and test functions
    % for the given norms.
    %
    % Parameters:
    %   usedAnsatzNorm: which norm to use for the computation of the
    %     normalization constants for the ansatz functions. @type string
    %   usedTestNorm: which norm to use for the computation of the normalization
    %     constants for the test functions. @type string
    % precomputeNormalizationConstants(obj, usedAnsatzNorm, usedTestNorm);

    % Normalization constant for a given ansatz function.
    %
    % Returns the normalization constant for an ansatz basis function given
    % through the indexes of the spatial and temporal basis functions. As
    % normalization is sometimes desired for different norms, the constants have
    % to be precomputed.
    %
    % See also:
    %   precomputeNormalizationConstants
    %
    % Parameters:
    %   jdx: index of spatial basis function @type integer
    %   kdx: index of temporal basis function @type integer
    %
    % Return values:
    %   val: normalization constant @type double
    % val = normalizationConstForAnsatzFunc(obj, jdx, kdx);

    val = normOfAnsatzFunc(obj, jdx, kdx);

    val = normOfTestFunc(obj, ldx, mdx, ndx);

    % Computes the normalization coefficients if normalization is used.
    precomputeNormalization();

    % Normalization constant for a given test function.
    %
    % Returns the normalization constant for an test basis function given
    % through the indexes of the spatial and temporal basis functions. As
    % normalization is sometimes desired for different norms, the constants have
    % to be precomputed.
    %
    % See also:
    %   precomputeNormalizationConstants
    %
    % Parameters:
    %   jdx: index of spatial basis function @type integer
    %   kdx: index of temporal basis function @type integer
    %   ndx: index of spatial basis function (initial condition) @type integer
    %
    % Return values:
    %   val: normalization constant @type double
    % val = normalizationConstForTestFunc(obj, ldx, mdx, ndx);
  end

  methods

    %% Custom getters and setter

    function dim = get.nAnsatzDim(obj)
      % Calculate the dimension of the ansatz space.

      dim = obj.nAnsatzSpatial * obj.nAnsatzTemporal;
    end

    function dim = get.nTestDim(obj)
      % Calculate the dimension of the test space.

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
      %   nTestSpatialIC: Number of sine basis functions for the initial
      %     condition @type integer

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

    function val = solutionFuncFromCoeffs(obj, solutionCoeffs, t, x, tspan)
      % Construct a function handle of the solution.
      %
      % Use it to define a solution function in (t, x) or evaluate the solution
      % directly for given t and x grids.
      %
      % Parameters:
      %   solutionCoeffs: solution vector of the linear system @type colvec
      %   t: temporal variable @type vector
      %   x: spatial variable @type vector
      %   tspan: custom temporal interval @type vector @default obj.tspan
      %
      % Return values:
      %   solfun: function handle of the solution function @type function_handle

      if nargin == 4
        tspan = obj.tspan;
      end

      val = zeros(size(t, 1), size(t, 2));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = cell(obj.nAnsatzSpatial, 1);
      for jdx = 1:obj.nAnsatzSpatial
        spatialValues{jdx} = obj.spatialBasisFunc(jdx, x);
      end
      temporalValues = cell(obj.nAnsatzTemporal, 1);
      for kdx = 1:obj.nAnsatzTemporal
        temporalValues{kdx} = obj.temporalBasisFunc(kdx, t, tspan);
      end

      for jdx = 1:obj.nAnsatzSpatial
        for kdx = 1:obj.nAnsatzTemporal
        % Get the right coefficient
        pos = (jdx - 1) * obj.nAnsatzTemporal + kdx;

        % @todo normalization constant
        if obj.useNormalization
          normAnsatz = obj.AnsatzNormDiag(pos, pos);
        else
          normAnsatz = 1;
        end

        % evaluate the corresponding basis functions
        val = val + solutionCoeffs(pos) * spatialValues{jdx} .* ...
          temporalValues{kdx} / normAnsatz;
        end
      end
    end

  end

end
