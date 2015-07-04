classdef SpaceTimeAssemblyAbstract < handle
  % Interface specification for the assembly of the space-time structures.
  %
  % The main purpose of this class is to provide a common core of methods which
  % are essential for the assembly of the needed matrices and vectors to solve
  % the space-time variational problem. The assembly of the space-time
  % structures relies on the tensor-product representation of the ansatz and
  % test functions which allows the computation of temporal and spatial matrices
  % and vectors.
  %
  % Further the assembly is divided into field independent and field dependent
  % parts which enables an easier use of this methods for the reduced basis
  % methods.

  properties
    % spatial interval @type vector
    xspan = [0 1];

    % temporal interval @type vector
    tspan = [0 1];
  end % properties

  methods(Abstract)
    % Evaluate the solution
    solval = evaluateSolution(obj, solvec, tgrid, xgrid);

    % Assemble the field independent part of the matrix
    assembleFieldIndependentMatrix(obj);

    % Assemble the load vector
    assembleLoadVector(obj);
  end % abstract methods

  methods

    function set.xspan(obj, val)
      % Setter for the spatial interval.
      %
      % Checks if the given value is a interval whose left end is null. This
      % condition is important for the chosen spatial basis functions.
      %
      % Parameters:
      %   val: given spatial interval @type rowvec

      if size(val, 1) ~= 1 || size(val, 2) ~= 2
        error('Given value is not a row vector of size 2!');
      elseif val(1) ~= 0
        error('The spatial interval must start at 0!');
      elseif val(2) <= val(1)
        error('The given spatial interval is empty or wrongly oriented!');
      end

      obj.xspan = val;
    end % set.xspan

  end % methods

end % classdef
