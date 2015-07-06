classdef SolverFourierLegendre < SolverAbstract
  % Solver for the propagators based on a Fourier base in space and Legendre
  % polynomials in time.
  %
  % @deprecated Not fully supported, because Legendre polynomials really aren't
  %     that suitable for the given task if there's more than one field.
  %
  % See also:
  %   SpatialAssemblyFourier TemporalAssemblyLegendre

  methods

    function obj = SolverFourierLegendre()
      % Constructor for this class

      % first things first: calls the superclass constructor, even if it doesn't
      % exist
      obj@SolverAbstract();

      % create the spatial and temporal assembly objects
      obj.spatial  = SpatialAssemblyFourier();
      obj.temporal = TemporalAssemblyLegendre();
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
      % Solve the forward propagator.
      %
      % Parameters:
      %   fieldCoefficients: cellarray of vectors that hold the coefficients for
      %     the field series expansions. @type cell
      %
      % Return values:
      %   solvec: coefficient vector of the solution in the trial space
      %     @type vector

      % first we set up the preconditioners.
      % attention: the inverse of these matrices will be multiplied with the
      % system matrix, not the matrices itself!
      % left side:
      % Pl = obj.TeNorm;
      Pl = speye(obj.nTestDim);
      % right side:
      % Pr = obj.TrNorm;
      Pr = speye(obj.nTrialDim);

      % sum the field dependent parts for the given series expansion
      % coefficients
      LhsFD = sparse(obj.nTestDim, obj.nTrialDim);
      % iterate over the fields
      for fdx = 1:obj.nFields
        % and now over the coefficients
        for cdx = 1:obj.nFieldCoeffs
          LhsFD = LhsFD + fieldCoefficients{fdx}(cdx) * obj.LhsFD{cdx, fdx};
        end
      end

      % apply the preconditioners to the system matrix
      LhsPre = (Pl \ (obj.LhsFI.F + LhsFD)) / Pr;

      % compute the right hand side load vector
      %| @todo implement a generalized case, this only allows uniform one initial
      %| condition and no source function.
      Rhs = sparse(obj.nTestDim, 1);
      Rhs(obj.nTestS * obj.nTestT + 1) = obj.xspan(2);

      % apply the left preconditioner
      RhsPre = Pl \ Rhs;

      % now solve the linear system
      solvecPre = LhsPre \ RhsPre;

      % and revert the preconditioning
      solvec = Pr \ solvecPre;
    end

    function solvec = solveBackward(obj, fieldCoefficients)
      % Solve the backward propagator.
      %
      % Parameters:
      %   fieldCoefficients: cellarray of vectors that hold the coefficients for
      %     the field series expansions. @type cell
      %
      % Return values:
      %   solvec: coefficient vector of the solution in the trial space
      %     @type vector

      % first we set up the preconditioners.
      % attention: the inverse of these matrices will be multiplied with the
      % system matrix, not the matrices itself!
      % left side:
      % Pl = obj.TeNorm;
      Pl = speye(obj.nTestDim);
      % right side:
      % Pr = obj.TrNorm;
      Pr = speye(obj.nTrialDim);

      % sum the field dependent parts for the given series expansion
      % coefficients
      LhsFD = sparse(obj.nTestDim, obj.nTrialDim);
      % iterate over the fields
      for fdx = 1:obj.nFields
        % and now over the coefficients
        for cdx = 1:obj.nFieldCoeffs
          LhsFD = LhsFD + fieldCoefficients{fdx}(cdx) * obj.LhsFD{cdx, fdx};
        end
      end

      % apply the preconditioners to the system matrix
      LhsPre = (Pl \ (obj.LhsFI.B + LhsFD)) / Pr;

      % compute the right hand side load vector
      %| @todo implement a generalized case, this only allows uniform one initial
      %| condition and no source function.
      Rhs = sparse(obj.nTestDim, 1);
      Rhs(obj.nTestS * obj.nTestT + 1) = obj.xspan(2);

      % apply the left preconditioner
      RhsPre = Pl \ Rhs;

      % now solve the linear system
      solvecPre = LhsPre \ RhsPre;

      % and revert the preconditioning
      solvec = Pr \ solvecPre;
    end

    function solval = evaluateSolution(obj, solvec, tgrid, xgrid)
      % Evaluate the solution.
      %
      % Parameters:
      %   solvec: coefficient vector of the solution in the trial space.
      %     @type vector.
      %   tgrid: temporal grid @type vector
      %   xgrid: spatial grid @type vector
      %
      % Return values:
      %   solval: values of the solution on the given grids @type matrix

      % first generate the meshgrid
      [tmesh, xmesh] = meshgrid(tgrid, xgrid);
      % and the matrix that will hold the solution values
      solval = zeros(size(tmesh, 1), size(tmesh, 2));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = cell(obj.nTrialS, 1);
      for jdx = 1:obj.nTrialS
        spatialValues{jdx} = obj.spatial.basisFunc(jdx, xmesh);
      end
      temporalValues = cell(obj.nTrialT, 1);
      for kdx = 1:obj.nTrialT
        temporalValues{kdx} = obj.temporal.basisFunc(kdx, tmesh);
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

  methods(Access = 'protected')

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
      % Assemble the field independent part of the space time system matrix.
      %
      % Warning:
      %   Only works for one field!
      %
      % Return values:
      %   FD: 2d cell array containing the part of the space time system matrix
      %     which corresponds to the (basis function, field) index pair
      %     @type cellarray
      %
      % @deprecated

      assert(obj.nFields == 1);

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

  end % protected methods

end % classdef
