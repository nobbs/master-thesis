classdef ProblemData < handle
  % Convenience class that stores the model variables.

  properties
    % Multiplicative constant of the Laplacian. @type double
    laplacian;
    % Added field offset. @type double
    offset

    % Spatial interval span. @type vector
    xspan;
    % Spatial interval grid. @type vector
    xgrid;
    % Temporal interval span. @type vector
    tspan;
    % Temporal interval grid. @type vector
    tgrid;

    % Points in time at which the field switches occur. @type vector
    f;
    % Number of basis functions to use per field. @type integer
    nC;
    % Index of basis functions to use. @type integer
    seriesIdx;
    % Number of fields @type integer
    nF;
  end

  properties(Dependent)
    % Total number of parameters. @type integer
    nP;
  end

  methods
    function val = get.nP(obj)
      % Total number of parameters. @type integer
      val = obj.nC * obj.nF;
    end
  end

end
