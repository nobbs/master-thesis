classdef SolverNodal < SolverAbstract
  % Solver for the propagators based on a Fourier base in space and piecewise
  % linear nodal functions in time.

  properties
    % Level of refinement for the temporal part of the test space. @type integer
    tref;

    % Toggles whether temporal grid refinement should be used for the test
    % space. Setting this to true guarantees stability! @type logical
    % useRefinement;
  end % public properties

  methods

    function obj = SolverNodal(pd, spatial, temporal)
      % Constructor for this class
      %
      % Parameters:
      %   pd: reference to the problem data object @type ProblemData
      %   spatial: spatial assembly object @type SpatialAssemblyAbstract
      %   temporal: nodal temporal assembly object @type TemporalAssemblyNodal

      assert(isa(temporal, 'TemporalAssemblyNodal'));

      % first things first: calls the superclass constructor, even if it doesn't
      % exist
      obj@SolverAbstract(pd, spatial, temporal);

      % call the prepare method to conclude the setup
      obj.prepare();
    end

    function prepare(obj)
      % Prepare the solver.
      %
      % This entails the assembly of the needed temporal and spatial structures
      % and the following combination via kronecker products to the space-time
      % structures.

      % compute all the needed space time structures
      obj.Lhs        = cell(obj.nQb, 1);
      obj.Lhs{1}     = obj.spacetimeStiffnessMatrix();
      obj.Lhs(2:end) = obj.spacetimeFieldDependentFourier();
      obj.TrNorm     = obj.spacetimeTrialNorm();
      obj.TeNorm     = obj.spacetimeTestNorm();
    end

    function solvec = solve(obj, param)
      % Solve the propagator.
      %
      % Parameters:
      %   param: matrix of the field series expansion coefficients. each column
      %     represents a field. @type matrix
      %
      % Return values:
      %   solvec: coefficient vector of the solution in the trial space
      %     @type vector

      if ~obj.temporal.useRefinement
        % first we set up the preconditioners.
        % attention: the inverse of these matrices will be multiplied with the
        % system matrix, not the matrices itself!
        % left side:
        % Pl = obj.TeNorm;
        Pl = speye(obj.nTestDim);
        % right side:
        % Pr = obj.TrNorm;
        Pr = speye(obj.nTrialDim);

        % assemble the system matrix for the given field
        Lhs = obj.spacetimeSystemMatrix(param);

        % apply the preconditioners to the system matrix
        LhsPre = (Pl \ Lhs) / Pr;

        % compute the right hand side load vector
        Rhs = obj.spacetimeLoadVector();

        % apply the left preconditioner
        RhsPre = Pl \ Rhs;

        % now solve the linear system
        solvecPre = LhsPre \ RhsPre;

        % and revert the preconditioning
        solvec = Pr \ solvecPre;
      else
        % BTN−1Bu = BTN−1b

        % first we set up the preconditioners.
        % attention: the inverse of these matrices will be multiplied with the
        % system matrix, not the matrices itself!
        % left side:
        % Pl = obj.TeNorm;
        % Pl = speye(obj.nTestDim);
        % right side:
        % Pr = obj.TrNorm;
        % Pr = speye(obj.nTrialDim);

        % assemble the system matrix for the given field
        B = obj.spacetimeSystemMatrix(param);

        Lhs = B.' * (obj.TeNorm \ B);

        % apply the preconditioners to the system matrix
        LhsPre = Lhs;
        % LhsPre = obj.TrNorm \ Lhs;

        % compute the right hand side load vector
        Rhs = obj.spacetimeLoadVector();
        Rhs = (B.' * (obj.TeNorm \ Rhs));
        % apply the left preconditioner
        RhsPre = Rhs;
        % RhsPre = obj.TrNorm \ Rhs;

        % now solve the linear system
        solvecPre = LhsPre \ RhsPre;

        % and revert the preconditioning
        solvec = solvecPre;
        % solvec = Pr \ solvecPre;
      end
    end

    function solval = evaluateSolution(obj, solvec)
      % Evaluate the solution.
      %
      % Parameters:
      %   solvec: coefficient vector of the solution in the trial space.
      %     @type vector.
      %
      % Return values:
      %   solval: values of the solution on the given grids @type matrix

      % create the matrix that will hold the solution values
      solval = zeros(length(obj.pd.tgrid), length(obj.pd.xgrid));

      % precompute the spatial and temporal basis functions for the given grids
      spatialValues = zeros(length(obj.pd.xgrid), obj.spatial.nTrial);
      for jdx = 1:obj.spatial.nTrial
        spatialValues(:, jdx) = obj.spatial.basisFunc(jdx, obj.pd.xgrid).';
      end

      % we evaluate the solution in two steps. the inner loop adds up all the
      % temporal evaluations that share the same spatial basis function as a
      % factor. the outer loop then multiplies this with the spatial component.
      for kdx = 1:obj.temporal.nTrial
        posStart       = (kdx - 1) * obj.spatial.nTrial + 1;
        posEnd         = kdx * obj.spatial.nTrial;
        solval(:, kdx) = spatialValues * solvec(posStart:posEnd);
      end

    end

  end

  methods%(Access = 'protected')
    function M = spacetimeStiffnessMatrix(obj)
      % Assemble the field independent part of the space time stiffness matrix.
      %
      % Return values:
      %   M: propagator field independent stiffness matrix @type matrix

      % temporal mass matrices
      MtAT   = obj.temporal.massMatrix('both');
      % temporal "half stiffness" matrix
      CtAT   = obj.temporal.halfStiffnessMatrix;

      % temporal propagation vectors. this and the sign in front of the time
      % derivative "half stiffness" matrix is the only difference between the
      % forward and the backward propagator!
      if obj.isForward
        et = obj.temporal.forwardInitVector;
        tdsign = 1;
      else
        et = obj.temporal.backwardInitVector;
        tdsign = -1;
      end

      % spatial mass matrices
      MxAT   = obj.spatial.massMatrix;
      MxATIc = obj.spatial.massMatrix;
      % spatial stiffness matrices
      AxAT   = obj.spatial.stiffnessMatrix;

      % compute the parts of the space time variational form as kronecker products
      timeDerivate = kron(CtAT, MxAT);
      laplacian    = kron(MtAT, AxAT);
      offset       = kron(MtAT, MxAT);
      initialCond  = kron(et, MxATIc);

      % set the upper block of the system matrix
      M = tdsign * timeDerivate + obj.pd.laplacian * laplacian + obj.pd.offset * offset;

      % resize it, so that we can add the lower block
      M(obj.nTestDim, obj.nTrialDim) = 0;

      % add the lower block which is responsible for the propagation of the
      % initial condition
      M((obj.temporal.nTest * obj.spatial.nTest + 1):end, :) = initialCond;
    end

    function FD = spacetimeFieldDependentFourier(obj)
      % Assemble the field independent part of the space time system matrix.
      %
      % Return values:
      %   FD: cell array containing the part of the space time system matrix
      %     which corresponds to the (basis func + field * num of basis func)
      %     index @type cellarray

      % spatial field dependent stuff
      FDx  = obj.spatial.fieldDependentFourier;

      % set up the needed cell array
      FD = cell(obj.pd.nP, 1);

      % get start and end points of the time interval parts
      tpoints = [obj.pd.tspan(1), obj.pd.f, obj.pd.tspan(2)];

      % iterate over the fields
      for fdx = 1:obj.pd.nF
        % get the time grid point indexes which are relevant for this field
        span = find(tpoints(fdx) <= obj.pd.tgrid & obj.pd.tgrid < tpoints(fdx + 1));

        % temporal mass matrices for the given interval part
        %| @todo add back the field support!
        MtAT = obj.temporal.massMatrix('both', [span(1), span(end)]);
        % MtAT = obj.temporal.massMatrix('both');

        % and iterate over the coefficients
        for cdx = 1:obj.pd.nC
          pos = (fdx - 1) * obj.pd.nC + cdx;
          % assemble
          FD{pos} = kron(MtAT, FDx{cdx});
          % and resize
          FD{pos}(obj.nTestDim, obj.nTrialDim) = 0;
        end
      end
    end

    function F = spacetimeLoadVector(obj)
      % Assemble the load vector for the space time system.
      % @todo generalize!

      F = sparse(obj.nTestDim, 1);
      F(obj.spatial.nTest * obj.temporal.nTest + 1) = obj.pd.xspan(2);
    end

    function M = spacetimeSystemMatrix(obj, param)
      % Assemble the system matrix for the given field coefficients.
      %
      % Parameters:
      %   param: vector of the field series expansion coefficients @type vector
      %
      % Return values:
      %   M: space time system matrix. @type matrix

      % sum the field dependent parts for the given series expansion
      % coefficients
      M = sparse(obj.nTestDim, obj.nTrialDim);

      % pad with one for the field independent part
      param = [1; shiftdim(param)];

      % add up the parts of the system matrix
      for pdx = 1:obj.nQb
        M = M + param(pdx) * obj.Lhs{pdx};
      end
    end

    function M = spacetimeTrialNorm(obj)
      % Assemble the matrix for the discrete norm on the trial space.
      %
      % Return values:
      %   M: gramian of the trial space norm @type matrix

      % temporal mass matrices
      MtAA = obj.temporal.massMatrix('trial');
      % temporal stiffness matrix
      AtAA = obj.temporal.stiffnessMatrix;
      % spatial mass matrices
      MxAA = obj.spatial.massMatrix;
      % spatial stiffness matrices
      AxAA = obj.spatial.stiffnessMatrix;

      % compute the norm for `L_2(I; V)` and `L_2(I; V')`
      M = kron(MtAA, MxAA + AxAA) + kron(AtAA, MxAA * ((MxAA + AxAA) \ MxAA));
    end

    function M = spacetimeTestNorm(obj)
      % Assemble the matrix for the discrete norm on the test space.
      %
      % Return values:
      %   M: gramian of the test space norm @type matrix

      % temporal mass matrices
      MtTT     = obj.temporal.massMatrix('test');
      % spatial mass matrices
      MxTT     = obj.spatial.massMatrix;
      % spatial stiffness matrices
      AxTT     = obj.spatial.stiffnessMatrix;

      % compute the two blocks of the norm
      A = kron(MtTT, MxTT + AxTT);
      B = MxTT;

      % compute the corresponding indices
      iA = 1:(obj.spatial.nTest * obj.temporal.nTest);
      iB = (obj.spatial.nTest * obj.temporal.nTest + 1):obj.nTestDim;

      % place the blocks
      M = sparse(obj.nTestDim, obj.nTestDim);
      M(iA, iA) = A;
      M(iB, iB) = B;
    end

  end % protected methods

end % classdef
