classdef ReducedBasisSolver < handle
  % Solver based on the reduced basis method.


  properties
    % underlying galerkin solver object
    solver;
  end

  methods
    function prepare(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setParameterData(obj)
      % Set parameter data.
      %
      % Used to set settings like number of parameters, ranges for the parameters
      % etc.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setParameterTrainingSet(obj)
      % Set parameter training set.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setDecisionAlgorithm(obj)
      % Set decision algorithm
      %
      % Used to choose the algorithm that will be used to determine which
      % parameter should be added next to the reduced basis.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function resetReducedBasis(obj)
      % Reset the reduced basis
      %
      % Clear the reduced basis. Won't delete the already assembled matrices and
      % vectors of galerkin method of the "truth" solution.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function train(obj)
      % Train reduced basis
      %
      % "Trains" the reduced basis method by evaluating the error estimator for a
      % lot of points of the training set and choosing a "worst" parameter which
      % will be added to the reduced basis.
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end
  end

  % Mathy stuff
  methods
    function estimateError(obj)
      % Error estimator
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function infSupLowerBound(obj)
      % Calculate lower bound for the inf-sup-constant
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function continuityUpperBound(obj)
      % Calculate upper bound for the continuity constant
      %
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    % more mathy stuff ...
  end

end
