classdef GalerkinSolver < handle
  % Galerkin based solver for the propagators.
  %
  % Theoretically this class should be able to handle different kinds of
  % underlying systems of basis functions for the ansatz and test subspace.
  % Theoretically. Right now it's only able to use the combination of Fourier
  % type spatial basis functions and Legendre polynomial temporal basis
  % functions.
  %
  % @todo Not fully implemented.

  properties
    % span of the spatial interval @type vector
    xspan            = [0, 1];
    % span of the temporal interval @type vector
    tspan            = [0, 1];
    % points in the time interval at which a switch of the field takes place @type vector
    breakpoints      = [];
    % multiplicative factor for the Laplacian @type double
    coeffLaplacian   = 1;
    % additive field-offset `\mu`. @type double
    coeffOffset      = 0;
    % Toggle if the ansatz and test functions should be normalized @type logical
    useNormalization = true;
    % number of coefficients of the field series expansions @type integer
    nFieldCoeffs     = 0;
  end % properties

  properties%(Access = 'private')
    % Cell array that holds the assembly objects and already assembled
    % structures. Each entry corresponds to a continuous part of the time
    % interval and an external field. @type cellarray
    parts;

    % number of basis functions for the underlying assembly classes
    nAnsatzSpatial;
    nAnsatzTemporal;
    useQuadraticSystem;
    nTestSpatial;
    nTestTemporal;
    nTestSpatialIC;
  end % properties private

  properties(Dependent)
    % total number of fields (or number of switches plus one)
    nFields;
  end % properties dependent

  methods

    % Constructor

    function obj = GalerkinSolver(varargin)
      % Constructor for the GalerkinSolver class.
      %
      % Create an object of the Galerkin Solver class with sensible default
      % values for not-so-important parameters. If you really want to use this,
      % then you obviously should set them to the values you have / need / want.
      %
      % Warning:
      %   If you want to set the number of the test subspace, then you've to set
      %   the three values nTestSpatial, nTestTemporal, nTestSpatialIC. If you
      %   don't set all three, then the set values are ignored and automatically
      %   chosen so that the resulting linear system is quadratic.
      %
      % Parameters:
      %   varargin: variable number of input parameters. These are in detail:
      %   breakpoints: the points in time at which the fields change @type
      %     vector
      %   Laplacian: scalar coefficient that is multiplied with the laplacian
      %     @type double
      %   nAnsatzSpatial: number of spatial basis functions for the ansatz
      %     subspace @type integer
      %   nAnsatzTemporal: number of temporal basis functions for the ansatz
      %     subspace @type integer
      %   nFieldCoeffs: number of spatial basis functions for the field series
      %     expansion @type integer
      %   nTestSpatial: number of spatial basis functions for the test subspace
      %     @type integer
      %   nTestTemporal: number of temporal basis functions for the test
      %     subspace @type integer
      %   nTestSpatialIC: number of spatial initial condition basis functions
      %     for the test subspace @type integer

      % set up the inputparser
      p = inputParser;

      % set up required parameters
      addRequired(p, 'breakpoints', @(x) isempty(x) || (isvector(x) & isnumeric(x)));
      addRequired(p, 'Laplacian', @(x) isnumeric(x) && isscalar(x));
      addRequired(p, 'nAnsatzSpatial', @(x) isnumeric(x) && isscalar(x));
      addRequired(p, 'nAnsatzTemporal', @(x) isnumeric(x) && isscalar(x));

      % set up optional parameters
      addOptional(p, 'nFieldCoeffs', 0, @(x) isnumeric(x) && isscalar(x));
      addOptional(p, 'nTestSpatial', 1, @(x) isnumeric(x) && isscalar(x));
      addOptional(p, 'nTestTemporal', 1, @(x) isnumeric(x) && isscalar(x));
      addOptional(p, 'nTestSpatialIC', 1, @(x) isnumeric(x) && isscalar(x));

      % set up optional optional parameters
      addParameter(p, 'useNormalization', obj.useNormalization, @islogical);
      addParameter(p, 'Offset', obj.coeffOffset, @(x) isnumeric(x) && isscalar(x));
      addParameter(p, 'tspan', obj.tspan, @(x) isnumeric(x) && isvector(x));
      addParameter(p, 'xspan', obj.xspan, @(x) isnumeric(x) && isvector(x));

      % parse input
      parse(p, varargin{:});

      % and now handle the input. first we forward the required inputs.
      obj.breakpoints     = p.Results.breakpoints;
      obj.coeffLaplacian  = p.Results.Laplacian;
      obj.nAnsatzSpatial  = p.Results.nAnsatzSpatial;
      obj.nAnsatzTemporal = p.Results.nAnsatzTemporal;

      % now we handle the optional test subspace size input.
      obj.nFieldCoeffs = p.Results.nFieldCoeffs;

      % if not all three optional arguments are set, we fall back to the default
      if any(strcmp('nTestSpatial', p.UsingDefaults) + ...
          strcmp('nTestTemporal', p.UsingDefaults) + ...
          strcmp('nTestSpatialIC', p.UsingDefaults)),
        obj.useQuadraticSystem = true;
      else
        obj.useQuadraticSystem = false;
        obj.nTestSpatial       = p.Results.nTestSpatial;
        obj.nTestTemporal      = p.Results.nTestTemporal;
        obj.nTestSpatialIC     = p.Results.nTestSpatialIC;
      end

      % and finally we'll handle the optional optional inputs
      obj.useNormalization = p.Results.useNormalization;
      obj.coeffOffset      = p.Results.Offset;
      obj.tspan            = p.Results.tspan;
      obj.xspan            = p.Results.xspan;
    end % GalerkinSolver

    function val = get.nFields(obj)
      val = length(obj.breakpoints) + 1;
    end % get.nFields


    %% Preparation

    function assemble(obj)
      % Creates all the structures needed to do anything.

      % as we want to solve the given pde for two fields, we will create two
      % assembly objects, one for each part of the temporal decomposition.
      parts = cell(obj.nFields, 1);

      % get the endpoints of the temporal decomposition
      tpoints = [obj.tspan(1), obj.breakpoints, obj.tspan(2)];

      % and now iterate over the parts of the decomposition and create the
      % needed assembly object and structures
      for fdx = 1:obj.nFields
        % create an object of the underlying assembly class
        parts{fdx}.assembly = AssemblyFourierLegendre();

        % set the number of used basis functions
        parts{fdx}.assembly.setNumberOfAnsatzFuncs(obj.nAnsatzSpatial, ...
          obj.nAnsatzTemporal);
        % if obj.useQuadraticSystem
          parts{fdx}.assembly.setNumberOfTestFuncsFromAnsatzFuncs();
        % else
        %   parts{fdx}.assembly.setNumberOfTestFuncs(obj.nTestSpatial, ...
        %     obj.nTestTemporal, obj.nTestSpatialIC);
        % end

        % set the spatial and temporal intervals
        parts{fdx}.assembly.xspan = obj.xspan;
        parts{fdx}.assembly.tspan = [tpoints(fdx), tpoints(fdx + 1)];

        % set whether we want to normalize the used ansatz and test basis
        % functions and precompute the needed norms for the normalization (has to
        % be called even if we don't want to normalize!)
        parts{fdx}.assembly.useNormalization = obj.useNormalization;
        parts{fdx}.assembly.precomputeNormalization();

        % assemble the field independent parts of the stiffness matrix
        parts{fdx}.matFI = parts{fdx}.assembly.assembleFieldIndependentMatrix();

        % and now assemble the field dependent part
        parts{fdx}.matFD = parts{fdx}.assembly.assembleFieldDependentMatrixForFourierSeries(obj.nFieldCoeffs);
      end

      % save the created stuff
      obj.parts = parts;
    end % assemble

    function solutionCoeffs = solveForward(obj, coefficients)
      % Solve the forward propagator for the given coefficients of series
      % expansions of the used external fields.
      %
      % Parameters:
      %   coefficients: cellarray of vectors with coefficients of the fields
      %
      % Return values:
      %   solutionCoeffs: coefficients of the solution expanded in the ansatz
      %     subspace basis functions @type vector

      solutionCoeffs = cell(obj.nFields, 1);

      % iterate over all the external fields
      for fdx = 1:obj.nFields
        matricesFI = obj.parts{fdx}.matFI;
        matricesFD = obj.parts{fdx}.matFD;

        % add up the field independent part of the stiffness Matrix
        LhsFI = matricesFI.TiD + obj.coeffLaplacian * matricesFI.Lap + ...
          obj.coeffOffset * matricesFI.Off + matricesFI.ICF;

        % add up the field dependent part of the stiffness matrix
        LhsFD = sparse(size(LhsFI, 1), size(LhsFI, 2));
        for cdx = 1:obj.nFieldCoeffs
          LhsFD = LhsFD + coefficients{fdx}(cdx) * matricesFD{cdx};
        end

        % construct the right hand side of the equation system. if this is the
        % first field, we start with the constant value of one as initial
        % conditions and the solution of the last interval part else.
        if fdx == 1
          Rhs = obj.parts{fdx}.assembly.assembleVectorOnes();
        else
          if obj.useNormalization
            initialCoeffs = obj.parts{fdx - 1}.assembly.AnsatzNormDiag \ ...
              solutionCoeffs{fdx - 1};
          else
            initialCoeffs = solutionCoeffs{fdx - 1};
          end
          Rhs = obj.parts{fdx}.assembly.assembleVectorFromSolutionCoeffs(initialCoeffs);
        end

        % solve the linear system
        solutionCoeffs{fdx} = (LhsFI + LhsFD) \ Rhs;
      end
    end % solveForward

    function solutionCoeffs = solveBackward(obj, coefficients)
      % Solve the backward propagator for the given coefficients of series
      % expansions of the used external fields.
      %
      % Parameters:
      %   coefficients: cellarray of vectors with coefficients of the fields
      %
      % Return values:
      %   solutionCoeffs: coefficients of the solution expanded in the ansatz
      %     subspace basis functions @type vector

      solutionCoeffs = cell(obj.nFields, 1);

      % iterate over all the external fields
      for fdx = obj.nFields:-1:1
        matricesFI = obj.parts{fdx}.matFI;
        matricesFD = obj.parts{fdx}.matFD;

        % add up the field independent part of the stiffness Matrix
        LhsFI = - matricesFI.TiD + obj.coeffLaplacian * matricesFI.Lap + ...
          obj.coeffOffset * matricesFI.Off + matricesFI.ICB;

        % add up the field dependent part of the stiffness matrix
        LhsFD = sparse(size(LhsFI, 1), size(LhsFI, 2));
        for cdx = 1:obj.nFieldCoeffs
          LhsFD = LhsFD + coefficients{fdx}(cdx) * matricesFD{cdx};
        end

        % construct the right hand side of the equation system. if this is the
        % first field, we start with the constant value of one as initial
        % conditions and the solution of the last interval part else.
        if fdx == obj.nFields
          Rhs = obj.parts{fdx}.assembly.assembleVectorOnes();
        else
          if obj.useNormalization
            initialCoeffs = obj.parts{fdx + 1}.assembly.AnsatzNormDiag \ ...
              solutionCoeffs{fdx + 1};
          else
            initialCoeffs = solutionCoeffs{fdx + 1};
          end
          Rhs = obj.parts{fdx}.assembly.assembleVectorFromSolutionCoeffs(initialCoeffs, true);
        end

        % solve the linear system
        solutionCoeffs{fdx} = (LhsFI + LhsFD) \ Rhs;
      end
    end % solveBackward

    function soleval = solCoeffsToSolFunc(obj, solutionCoeffs, t, x)
      % Evaluate the solution function for the given solution coefficients and
      % the given temporal and spatial values.
      %
      % Parameters:
      %   solutionCoeffs: cellarray of coefficients of the solution
      %   t: points in the time interval @type vector
      %   x: points in the space interval @type vector
      %
      % Return values:
      %   solfun: function handle of the solution function @type function_handle
      %
      % @todo somehow we have to assure that the grid has always the same
      %     structure / orientation

      soleval = zeros(size(t, 1), size(t, 2));
      tpoints = [obj.tspan(1), obj.breakpoints, obj.tspan(2)];

      % get the corresponding part of the grids
      Idx = tpoints(1) <= t(1, :) & t(1, :) <= tpoints(2);
      tp = t(:, Idx);
      xp = x(:, Idx);
      % add the first field manually to get the initial conditions right
      soleval(:, Idx) = obj.parts{1}.assembly.solutionFuncFromCoeffs(solutionCoeffs{1}, tp, xp);

      % iterate over the parts of the temporal decomposition and add up the
      % corresponding solution evaluated at the given points
      for fdx = 2:obj.nFields
        % get the corresponding part of the grids
        Idx = tpoints(fdx) < t(1, :) & t(1, :) <= tpoints(fdx + 1);
        tp = t(:, Idx);
        xp = x(:, Idx);
        soleval(:, Idx) = obj.parts{fdx}.assembly.solutionFuncFromCoeffs(...
            solutionCoeffs{fdx}, tp, xp);
      end
    end % solCoeffsToSolFunc

  end % methods

end
