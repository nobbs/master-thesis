classdef SolverAbstract < handle
  % Common interface for the solver classes.
  % @deprecated

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
  end % properties

  properties(Dependent)
    % total number of fields (or number of switches plus one) @type integer
    nFields;

    % dimension of the trial space @type integer
    nTrialDim;
    % dimension of the test space @type integer
    nTestDim;
  end % dependent properties

  methods(Abstract)
    % Prepare the solver
    prepare(obj);

    % Solve the forward propagator.
    solveForward(obj, fieldCoefficients);

    % Solve the backward propagator.
    solveBackward(obj, fieldCoefficients);

    % Evaluate the solution
    evaluateSolution(obj, solvec);
  end % abstract methods

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
end % classdef
