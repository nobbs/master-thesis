classdef SolverFourierNodal < SolverAbstract
  % Solver for the propagators based on a Fourier base in space and piecewise
  % linear nodal functions in time.
  %
  % See also:
  %   SpatialAssemblyFourier TemporalAssemblyLinearConstant

  properties
    % Level of refinement for the temporal part of the test space. @type integer
    tref;

    % Toggles whether temporal grid refinement should be used for the test
    % space. Setting this to true guarantees stability! @type logical
    useRefinement;
  end % public properties

  methods

    function obj = SolverFourierNodal(pd, nspatial, useRefinement)
      % Constructor for this class
      %
      % Parameters:
      %   pd: problem data object @type ProblemData
      %   nspatial: number of spatial functions to use for trial and test space
      %     @type integer
      %   useRefinement: toggle whether refinement for the test space should be
      %     used. @type logical @default true

      % set default options
      if nargin == 2
        useRefinement = true;
      end

      % first things first: calls the superclass constructor, even if it doesn't
      % exist
      obj@SolverAbstract();

      % save the problem data object
      obj.pd = pd;

      % refinement settings
      obj.useRefinement = useRefinement;

      % forward the number of spatial functions
      obj.nTrialS  = nspatial;
      obj.nTestS   = nspatial;
      obj.nTestSic = nspatial;

      % create the spatial and temporal assembly objects
      obj.spatial  = SpatialAssemblyFourier();
      obj.temporal = TemporalAssemblyLinearConstant();

      % call the prepare method to conclude the setup
      obj.prepare();
    end

    function rhs = rhs(obj)
      % Calculate the load vector.
      % @todo generalize and move it!

      Rhs = sparse(obj.nTestDim, 1);
      Rhs(obj.nTestS * obj.nTestT + 1) = obj.pd.xspan(2);

      % apply the left preconditioner
      % RhsPre = Pl \ Rhs;
      rhs = Rhs;
    end

    function prepare(obj)
      % Prepare the solver.
      %
      % This entails the assembly of the needed temporal and spatial structures
      % and the following combination via kronecker products to the space-time
      % structures.

      % copy stuff from the problem data object
      obj.nTrialT = length(obj.pd.tgrid);
      obj.nTestT = (1 + (obj.useRefinement)) *(length(obj.pd.tgrid) - 1);

      % forward the must-have values to the assembly objects
      obj.spatial.xwidth = obj.pd.xspan(2) - obj.pd.xspan(1);
      obj.temporal.tgrid = obj.pd.tgrid;

      if obj.useRefinement
        obj.nTestT = 2 * (obj.nTrialT - 1);
        obj.tref = 1;
      else
        obj.tref = 0;
      end

      % and now compute all the needed space time structures
      obj.Lhs = cell(obj.nQb, 1);
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

      if ~obj.useRefinement
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
      spatialValues = zeros(length(obj.pd.xgrid), obj.nTrialS);
      for jdx = 1:obj.nTrialS
        spatialValues(:, jdx) = obj.spatial.basisFunc(jdx, obj.pd.xgrid).';
      end

      % we evaluate the solution in two steps. the inner loop adds up all the
      % temporal evaluations that share the same spatial basis function as a
      % factor. the outer loop then multiplies this with the spatial component.
      for kdx = 1:obj.nTrialT
        posStart = (kdx - 1) * obj.nTrialS + 1;
        posEnd = kdx * obj.nTrialS;
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
      MtAT   = obj.temporal.massMatrixBoth(obj.pd.tgrid, obj.tref);
      % temporal "half stiffness" matrix
      CtAT   = obj.temporal.halfStiffnessMatrix(obj.pd.tgrid, obj.tref);

      % temporal propagation vectors. this and the sign in front of the time
      % derivative "half stiffness" matrix is the only difference between the
      % forward and the backward propagator!
      if obj.isForward
        et = obj.temporal.forwardInitVector();
        tdsign = 1;
      else
        et = obj.temporal.backwardInitVector();
        tdsign = -1;
      end

      % spatial mass matrices
      MxAT   = obj.spatial.massMatrix(obj.nTrialS, obj.nTestS);
      MxATIc = obj.spatial.massMatrix(obj.nTrialS, obj.nTestSic);
      % spatial stiffness matrices
      AxAT   = obj.spatial.stiffnessMatrix(obj.nTrialS, obj.nTestS);

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
      M((obj.nTestT * obj.nTestS + 1):end, :) = initialCond;
    end

    function FD = spacetimeFieldDependentFourier(obj)
      % Assemble the field independent part of the space time system matrix.
      %
      % Return values:
      %   FD: cell array containing the part of the space time system matrix
      %     which corresponds to the (basis func + field * num of basis func)
      %     index @type cellarray

      % spatial field dependent stuff
      FDx  = obj.spatial.fieldDependentFourier(obj.nTrialS, obj.nTestS, obj.pd.nC);

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
        % MtAT = obj.temporal.massMatrix('both', obj.useRefinement, [span(1), span(end)]);
        MtAT = obj.temporal.massMatrixBoth(obj.pd.tgrid, obj.tref);

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
      F(obj.nTestS * obj.nTestT + 1) = obj.pd.xspan(2);
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
      MtAA = obj.temporal.massMatrixTrial();
      % temporal stiffness matrix
      AtAA = obj.temporal.stiffnessMatrix();
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
      MtTT     = obj.temporal.massMatrixTest(obj.pd.tgrid, obj.tref);
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
