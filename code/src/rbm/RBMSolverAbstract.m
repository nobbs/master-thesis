classdef RBMSolverAbstract < handle
  % Common interface for the reduced basis method solver classes.
  %
  % Provides common properties and methods so we don't have to define them for
  % every possible solver.

  properties
    % Number of basis functions to use for the field series expansion.
    % @type integer
    nC;
    % Points in time at which the field switches occur. @type vector
    breakpoints;

    trialshots;
    testshots;
  end % properties

  properties(Dependent)
    % Total number of fields. @type integer
    nFields;
  end % dependent properties

  properties%(Access = 'protected')
    % Underlying galerkin solver object @type SolverAbstract
    solver;
  end % protected properties

  methods(Abstract)
    % Prepare the RBM Solver.
    prepare(obj);

    % Offline phase.
    offlinePhase(obj);

    % Solve online for a given Parameter.
    % onlineSolve(obj);

    % Solve offline for a given Parameter.
    offlineSolve(obj);

    % Train the solver.
    % train(obj);

    % Delete the data gained through training.
    % resetTraining(obj);

    % Extend the RBM test space.
    % extendTestSpace(obj);

    % Estimate the current RBM-Truth-Error.
    % errorEstimate(obj);

    % Compute the Y-Norm of the residual.
    residual(obj);

    % Compute a lower bound for the inf-sup-constant.
    % infSupLowerBound(obj);

    % Compute a upper bound for the continuity constant.
    % continuityUpperBound(obj);

    % Evaluate a solution of the reduced basis method.
    evaluateSolutionRb(obj, solvec);

    % Evaluate a solution of the underlying galerkin solver.
    evaluateSolutionTruth(obj, solvec);
  end % abstract methods

  methods(Abstract)%, Access = 'protected')
    % Calculate the discrete inf-sup-constant of the rb system for the given
    % parameter.
    [beta, gamma] = calcDiscreteInfSupAndContinuityRB(obj, param);

    % Calculate the discrete inf-sup-constant of the truth system for the
    % given parameter.
    [beta, gamma] = calcDiscreteInfSupAndContinuityTruth(obj, param);
  end % private abstract methods

  methods % for dependent properties
    function val = get.nFields(obj)
      % Total number of fields.
      %
      % Return values:
      %   val: total number of fields @type integer

      val = length(obj.breakpoints) + 1;
    end
  end % methods for dependent properties

end % classdef
