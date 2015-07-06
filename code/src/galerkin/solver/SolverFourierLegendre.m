classdef SolverFourierLegendre < handle
  % @deprecated not fully supported!

  properties
    % number of spatial basis functions for the trial space @type integer
    nTrialS;
    % number of temporal basis functions for the trial space @type integer
    nTrialT;
    % number of spatial basis functions for the test space @type integer
    nTestS;
    % number of temporal basis functions for the test space @type integer
    nTestT;
    % number of spatial basis functions for the initial condition in the test
    % space @type integer
    nTestSic;

    % span of the spatial interval @type vector
    xspan;
    % span of the temporal interval @type vector
    tspan;

    % multiplicative factor for the Laplacian @type double
    cLaplacian;
    % additive field-offset `\mu`. @type double
    cOffset;

    % points in time at which the the field switch occurs @type vector
    breakpoints;
    % number of coefficients of the field series expansions @type integer
    nFieldCoeffs;

    % field dependent matrices in a cell array @type struct
    FDx;
  end

  properties(Dependent)
    % total number of fields (or number of switches plus one) @type integer
    nFields;

    % dimension of the trial space @type integer
    nTrialDim;
    % dimension of the test space @type integer
    nTestDim;
  end % dependent properties

  methods
    function val = get.nFields(obj)
      % Return the value of fields.
      val = length(obj.breakpoints) + 1;
    end

    function val = get.nTrialDim(obj)
      % Return the dimension of the trial space.
      val = obj.nTrialT * obj.nTrialS;
    end

    function val = get.nTestDim(obj)
      % Return the dimension of the test space.
      val = obj.nTestT * obj.nTestS + obj.nTestSic;
    end
  end

  properties(Access = 'private')
    % spatial assembly object @type SpatialAssemblyFourier
    spatial;

    % temporal assembly object @type TemporalAssemblyLegendre
    temporal;

    % struct that holds the field independent parts of the space time stiffness
    % matrix for the forward and backward propagator
    LhsFI;

    LhsFD;

    % trial norm matrix @type matrix
    TrNorm;

    % test norm matrix @type matrix
    TeNorm;
  end

  methods

    function obj = SolverFourierLegendre()
      % Constructor for this class

      % create the spatial and temporal assembly objects
      obj.spatial  = SpatialAssemblyFourier();
      obj.temporal = TemporalAssemblyLegendre();

      % create the needed structures
      obj.FDx    = struct();
    end

    function prepare(obj)
      % Prepare the solver.
      %
      % This entails the assembly of the needed temporal and spatial structures
      % and the following combination via kronecker products to the space-time
      % structures.

      % forward the current spatial and temporal intervals to the assembly
      % objects
      obj.spatial.xwidth = obj.xspan(2) - obj.xspan(1);
      obj.temporal.twidth = obj.tspan(2) - obj.tspan(1);

      % compute the space time structures;
      obj.LhsFI  = obj.spacetimeStiffnessMatrix();
      obj.LhsFD  = obj.spacetimeFieldDependentFourier();
      obj.TrNorm = obj.spacetimeTrialNorm();
      obj.TeNorm = obj.spacetimeTestNorm();
    end

    function solvec = solveForward(obj, fieldCoefficients)
      % @todo add field support


      % left side preconditioner (the inverse will be used)
      PL = sqrt(obj.TeNorm);
      % right side preconditioner (the inverse will be used)
      PR = sqrt(obj.TrNorm);
      % PR = speye(obj.nTrialDim);

      % sum the field dependent part

      LhsFD = sparse(obj.nTestDim, obj.nTrialDim);
      size(obj.LhsFD)
      for fdx = 1:obj.nFields
        for cdx = 1:obj.nFieldCoeffs
          LhsFD = LhsFD + fieldCoefficients{fdx}(cdx) * obj.LhsFD{cdx, fdx};
        end
      end

      % precondition the left hand side
      LhsPre = (PL \ (obj.LhsFI.F + LhsFD)) / PR;

      % compute the right hand side load vector
      Rhs = sparse(obj.nTestDim, 1);
      Rhs(obj.nTestS * obj.nTestT + 1) = obj.xspan(2)

      % cond(full(LhsPre))
      % return;

      % and consequently also precondition it
      RhsPre = PL \ Rhs;

      % solve the system
      solvecPre = LhsPre \ RhsPre;

      % and now revert the preconditioning
      solvec = PR \ solvecPre;
    end

    function solvec = solveBackward(obj, fieldCoefficients)
      % @todo add field support

      % left side preconditioner (the inverse will be used)
      PL = sqrt(obj.TeNorm);
      % right side preconditioner (the inverse will be used)
      % PR = obj.TrNorm;
      PR = speye(obj.nTrialDim);

      % precondition the left hand side
      LhsPre = (PL \ obj.LhsFI.B) / PR;

      % compute the right hand side load vector
      Rhs = sparse(obj.nTestDim, 1);
      Rhs(obj.nTestS * obj.nTestT + 1) = obj.xspan(2);

      % and consequently also precondition it
      RhsPre = PL \ Rhs;

      % solve the system
      solvecPre = LhsPre \ RhsPre;

      % and now revert the preconditioning
      solvec = PR \ solvecPre;
    end

    function solval = evaluateSolution(obj, solvec, tmesh, xmesh)

      solval = zeros(size(tmesh, 1), size(tmesh, 2));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = cell(obj.nTrialS, 1);
      for jdx = 1:obj.nTrialS
        spatialValues{jdx} = obj.spatialBasisFunc(jdx, xmesh);
      end
      temporalValues = cell(obj.nTrialT, 1);
      for kdx = 1:obj.nTrialT
        temporalValues{kdx} = obj.temporalBasisFunc(kdx, tmesh);
      end

      % we evaluate the solution in two steps. the inner loop adds up all the
      % temporal evaluations that share the same spatial basis function as a
      % factor. the outer loop then multiplies this with the spatial component.
      for jdx = 1:obj.nTrialS
        % temporalVal = zeros(size(tmesh, 1), size(tmesh, 2));
        for kdx = 1:obj.nTrialT
          % Get the right coefficient
          pos = (kdx - 1) * obj.nTrialS + jdx;

          % evaluate the corresponding basis functions
          solval = solval + solvec(pos) * temporalValues{kdx} .* spatialValues{jdx};
        end
        % solval = solval + temporalVal .* spatialValues{jdx};
      end

    end

  end

  methods(Access = 'private')

    function M = spacetimeStiffnessMatrix(obj)
      % Assemble the field independent part of the space time stiffness matrix.
      %
      % Return values:
      %   M: @type struct. In Detail:
      %      F: forward propagator field independent stiffness matrix
      %      B: backward propagator field independent stiffness matrix

      % temporal mass matrices
      MtAT   = obj.temporal.massMatrix(obj.nTrialT, obj.nTestT);
      % temporal "half stiffness" matrix
      CtAT   = obj.temporal.halfStiffnessMatrix(obj.nTrialT, obj.nTestT);
      % temporal propagation vectors
      etF    = obj.temporal.forwardInitVector(obj.nTrialT);
      etB    = obj.temporal.backwardInitVector(obj.nTrialT);
      % spatial mass matrices
      MxAT   = obj.spatial.massMatrix(obj.nTrialS, obj.nTestS);
      MxATIc = obj.spatial.massMatrix(obj.nTrialS, obj.nTestSic);
      % spatial stiffness matrices
      AxAT   = obj.spatial.stiffnessMatrix(obj.nTrialS, obj.nTestS);

      % compute the parts of the space time variational form as kronecker products
      timeDerivate = kron(CtAT, MxAT);
      laplacian    = kron(MtAT, AxAT);
      offset       = kron(MtAT, MxAT);
      initialF     = kron(etF, MxATIc);
      initialB     = kron(etB, MxATIc);

      M = struct();

      % set the upper block of the forward and backward propagator stiffness
      % matrices
      M.F =  timeDerivate + obj.cLaplacian * laplacian + obj.cOffset * offset;
      M.B = -timeDerivate + obj.cLaplacian * laplacian + obj.cOffset * offset;

      % resize both, so that we can add the lower block
      M.F(obj.nTestDim, obj.nTrialDim) = 0;
      M.B(obj.nTestDim, obj.nTrialDim) = 0;

      % add the lower block which is responsible for the propagation of the
      % initial condition
      M.F((obj.nTestT * obj.nTestS + 1):end, :) = initialF;
      M.B((obj.nTestT * obj.nTestS + 1):end, :) = initialB;
    end

    function FD = spacetimeFieldDependentFourier(obj)
      % @todo find out, how far the temporal matrix has to go...

      % temporal mass matrices
      MtAT = obj.temporal.massMatrix(obj.nTrialT, obj.nTestT);

      % spatial field dependent stuff
      FDx  = obj.spatial.fieldDependentFourier(obj.nTrialS, obj.nTestS, obj.nFieldCoeffs);

      % set up the needed cell array
      FD = cell(obj.nFieldCoeffs, obj.nFields);

      for fdx = 1:obj.nFields
        for cdx = 1:obj.nFieldCoeffs
          % assemble
          FD{cdx, fdx} = kron(MtAT, FDx{cdx});
          % resize
          FD{cdx, fdx}(obj.nTestDim, obj.nTrialDim) = 0;
        end
      end
    end

    function M = spacetimeTrialNorm(obj)
      % Assemble the matrix for the discrete norm on the trial space.
      %
      % Return values:
      %   M: gramian of the trial space norm @type matrix

      % temporal mass matrices
      MtAA = obj.temporal.massMatrix(obj.nTrialT, obj.nTrialT);
      % temporal stiffness matrix
      AtAA = obj.temporal.stiffnessMatrix(obj.nTrialT);
      % spatial mass matrices
      MxAA = obj.spatial.massMatrix(obj.nTrialS, obj.nTrialS);
      % spatial stiffness matrices
      AxAA = obj.spatial.stiffnessMatrix(obj.nTrialS, obj.nTrialS);

      % compute the norm for `L_2(I; V)` and `L_2(I; V')`
      M = kron(MtAA, MxAA + AxAA) + kron(AtAA, MxAA * ((MxAA + AxAA) \ MxAA));
    end

    function M = spacetimeTestNorm(obj)
      % Assemble the matrix for the discrete norm on the test space.
      %
      % Return values:
      %   M: gramian of the test space norm @type matrix

      % temporal mass matrices
      MtTT     = obj.temporal.massMatrix(obj.nTestT, obj.nTestT);
      % spatial mass matrices
      MxTT     = obj.spatial.massMatrix(obj.nTestS, obj.nTestS);
      MxTicTic = obj.spatial.massMatrix(obj.nTestSic, obj.nTestSic);
      % spatial stiffness matrices
      AxTT     = obj.spatial.stiffnessMatrix(obj.nTestS, obj.nTestS);

      % compute the two blocks of the norm
      A = kron(MtTT, MxTT + AxTT);
      B = MxTicTic;

      % compute the corresponding indices
      iA = 1:(obj.nTestS * obj.nTestT);
      iB = (obj.nTestS * obj.nTestT + 1):obj.nTestDim;

      % place the blocks
      M = sparse(obj.nTestDim, obj.nTestDim);
      M(iA, iA) = A;
      M(iB, iB) = B;
    end

    function val = spatialBasisFunc(obj, index, x)
      % Spatial basis functions.
      %
      % Evaluates the spatial basis function, Fourier functions, with the given
      % index for the given values of x. Can be used to define function handles
      % and numerical integration.
      %
      % In this case the spatial basis functions are of the type
      % ``\begin{cases}
      %   1, & i = 1\\
      %   \sin(\pi i x / L), & i~\text{even}\\
      %   \cos(\pi i x / L), & i~\text{odd and}~i > 1
      % \end{cases},`` where `i` corresponds to index and `L` is the width of
      % the spatial interval.
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   x: values in which the function should be evaluated @type matrix
      %
      % Return values:
      %   val: values of the basis function in x @type matrix

      if index == 1
        val = ones(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = sin(pi * index * x / obj.xspan(2));
      else
        val = cos(pi * (index - 1) * x / obj.xspan(2));
      end
    end

    function val = temporalBasisFunc(obj, index, t)
      % Temporal basis functions.
      %
      % Evaluates the temporal basis function, shifted Legendre polynomials,
      % with the given index for the given values of t. Can be used to define
      % function handles and numerical integration.
      %
      % Warning:
      %   The enumeration of index starts at 1, that means you'll get the
      %   Legendre polynomial of degree (index - 1)!
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   t: values in which the function should be evaluated @type matrix
      %   tspan: custom temporal interval @type vector @default obj.tspan
      %
      % Return values:
      %   val: values of the basis function in t @type matrix

      val = legendrePolynomial(t, index - 1, obj.tspan);
    end % spatialBasisFuncDerivative

  end % private methods

end % classdef
