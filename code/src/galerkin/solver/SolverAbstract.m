classdef SolverAbstract < handle
  % Common interface for the solver classes.
  %
  % Provides common properties and methods so we don't have to define them for
  % every possible solver.

  properties
    % Number of trial space spatial basis functions. @type integer
    nTrialS;
    % Number of trial space temporal basis functions. @type integer
    nTrialT;
    % Number of test space spatial basis functions. @type integer
    nTestS;
    % Number of test space temporal basis functions. @type integer
    nTestT;
    % Number of test space initial condition spatial basis functions.
    % @type integer
    nTestSic;

    % Span of the spatial interval. @type vector
    xspan;
    % Span of the temporal interval. @type vector
    tspan;

    % Direction of the propagator. @type logical
    isForward = true;

    % Multiplicative constant of the Laplace-Operator. @type double
    cLaplacian;
    % Added field offset `\mu`. @type double
    cOffset;

    % Points in time at which the field switches occur. @type vector
    breakpoints;
    % Number of basis functions to use for the field series expansion.
    % @type integer
    nC;

    % Total number of parameters @type integer
    nP;

    % Parts of the space time system matrix. The first entry is the field
    % independent part followed by the field dependent parts. @type cellarray
    Lhs;

    % As the right hand side of our system is the same every time, let's save
    % it. @type vector
    % @deprecated
    Rhs;

    % Norm of the trial space. @type matrix
    TrNorm;
    % Norm of the test space. @type matrix
    TeNorm;
  end % properties

  properties%(Access = 'protected')
    % Object responsible for the assembly of spatial structures.
    % @type SpatialAssemblyAbstract
    spatial;
    % Object responsible for the assembly of temporal structures.
    % @type TemporalAssemblyAbstract
    temporal;
  end % protected properties

  properties(Dependent)
    % Dimension of the trial space. @type integer
    nTrialDim;
    % Dimension of the test space. @type integer
    nTestDim;

    % Total number of fields. @type integer
    nFields;
    % Total number of bilinear forms @type integer
    nQb;
  end % dependent properties

  methods(Abstract)
    % Prepare the solver.
    %
    % This mainly consist of assembling the space time system matrix and the
    % matrices of the norms for the trial and test space.
    prepare(obj);

    % Solve the propagator.
    %
    % Parameters:
    %   param: matrix of the field series expansion coefficients. each column
    %     represents a field. @type matrix
    solve(obj, param);

    % Evaluate the solution.
    %
    % Parameters:
    %   solvec: coefficient vector of the solution in the trial space.
    %     @type vector.
    evaluateSolution(obj, solvec);
  end % abstract methods

  methods(Abstract)%, Access = 'protected')
    % Assemble the field independent part of the space time stiffness matrix.
    spacetimeStiffnessMatrix(obj);

    % Assemble the field independent part of the space time system matrix.
    spacetimeFieldDependentFourier(obj);

    % Assemble the matrix for the discrete norm on the trial space.
    spacetimeTrialNorm(obj);

    % Assemble the matrix for the discrete norm on the test space.
    spacetimeTestNorm(obj);
  end % abstract protected methods

  methods % for dependent properties
    function val = get.nFields(obj)
      % Total number of fields.
      %
      % Return values:
      %   val: total number of fields @type integer

      val = length(obj.breakpoints) + 1;
    end

    function val = get.nTrialDim(obj)
      % Dimension of the trial space.
      %
      % Return values:
      %   val: dimension of the trial space @type integer

      val = obj.nTrialT * obj.nTrialS;
    end

    function val = get.nTestDim(obj)
      % Dimension of the test space.
      %
      % Return values:
      %   val: dimension of the test space @type integer

      val = obj.nTestT * obj.nTestS + obj.nTestSic;
    end

    function val = get.nQb(obj)
      val = 1 + obj.nFields * obj.nC;
    end
  end % methods for dependent properties

  % methods(Access = 'protected')

  % end % protected methods

end % classdef
