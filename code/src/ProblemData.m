classdef ProblemData < handle
  % Convenience class that stores the model variables.

  properties
    % Multiplicative constant of the Laplacian. @type double
    laplacian;
    % Added field offset. @type double
    offset;

    % Initial condition given by a double value or a function handle.
    icfun;
    % Source term given by a double value, a spatial only function handle or a
    % space time function.
    sourcefun;

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
    % Toggle between sine and fourier series expansion for the fields (defaults
    % to fourier). @type logical
    useSineExpansion = false;
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
