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

    % Multiplicative constant of the Laplace-Operator. @type double
    cLaplacian;
    % Added field offset `\mu`. @type double
    cOffset;

    % Points in time at which the field switches occur. @type vector
    breakpoints;
    % Number of basis functions to use for the field series expansion.
    % @type integer
    nFieldCoeffs;
  end % properties

  properties(Access = 'protected')
    % Object responsible for the assembly of spatial structures.
    % @type SpatialAssemblyAbstract
    spatial;
    % Object responsible for the assembly of temporal structures.
    % @type TemporalAssemblyAbstract
    temporal;

    % Field independent parts of the space time system matrix. @type struct
    LhsFI;
    % Field dependent parts of the space time system matrix. @type cellarray
    LhsFD;

    % Norm of the trial space. @type matrix
    TrNorm;
    % Norm of the test space. @type matrix
    TeNorm;
  end % protected properties

  properties(Dependent)
    % Dimension of the trial space. @type integer
    nTrialDim;
    % Dimension of the test space. @type integer
    nTestDim;

    % Total number of fields. @type integer
    nFields;
  end % dependent properties

  methods(Abstract)
    % Prepare the solver.
    %
    % This mainly consist of assembling the space time system matrix and the
    % matrices of the norms for the trial and test space.
    prepare(obj);

    % Solve the forward propagator.
    %
    % Parameters:
    %   fieldCoefficients: cellarray of vectors that hold the coefficients for
    %     the field series expansions. @type cell
    solveForward(obj, fieldCoefficients);

    % Solve the backward propagator.
    %
    % Parameters:
    %   fieldCoefficients: cellarray of vectors that hold the coefficients for
    %     the field series expansions. @type cell
    solveBackward(obj, fieldCoefficients);

    % Evaluate the solution.
    %
    % Parameters:
    %   solvec: coefficient vector of the solution in the trial space.
    %     @type vector.
    evaluateSolution(obj, solvec);
  end % abstract methods

  methods(Abstract, Access = 'protected')
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
  end % methods for dependent properties
end % classdef
