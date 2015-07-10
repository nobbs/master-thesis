classdef RBMSolverAbstract < handle
  % Common interface for the reduced basis method solver classes.
  %
  % Provides common properties and methods so we don't have to define them for
  % every possible solver.

  properties
    % Reference to the given problem data object. @type ProblemData
    pd;
    % Underlying galerkin "truth" solver object @type SolverAbstract
    solver;
  end % properties

  properties(Access = 'protected')
    % Holds the truth solution vectors which form the basis of the reduced basis
    % trial space. @type matrix
    trialSnapshots;
  end % protected properties

  methods
    function obj = RBMSolverAbstract(problem)
      % Default constructor.
      %
      % Parameters:
      %   problem: reference to a problem data object. @type ProblemData

      % save the reference, nothing more
      obj.pd = problem;
    end
  end

  methods(Abstract)
    % Prepare the RBM Solver.
    prepare(obj);

    % Offline stage.
    offlineStage(obj);

    % Solve online for a given Parameter with the reduced basis solver.
    onlineSolve(obj, param);

    % Delete the data gained through training.
    resetTraining(obj);

    % Estimate the current RBM-Truth-Error.
    estimateError(obj);

    % Compute the Y-Norm of the residual.
    residual(obj);

    % Evaluate a solution of the reduced basis solver.
    evaluateSolutionRb(obj, solvec);

    % Evaluate a solution of the truth solver.
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

end % classdef
