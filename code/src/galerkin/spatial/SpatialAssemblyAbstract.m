classdef SpatialAssemblyAbstract < handle
  % Interface for the assembly of spatial structures.
  %
  % Needed for the construction of the space-time-structures, as these are
  % computed through Kronecker-products of temporal and spatial matrices and
  % vectors.
  %
  % See also:
  %   TemporalAssemblyAbstract

  properties
    % Reference to the problem data object @type ProblemData
    pd;
    % Number of basis functions in the trial space @type integer
    nTrial;
    % Number of basis functions in the test space @type integer
    nTest;
    % Number of basis functions in the test space for the initial condition
    % @type integer
    nTestIC;
  end

  methods
    function obj = SpatialAssemblyAbstract(pd, nTrial, nTest, nTestIC)
      % Default constructor
      %
      % Parameters:
      %   pd: Reference to the problem data object @type ProblemData
      %   nTrial: Number of basis functions in the trial space @type integer
      %   nTest: Number of basis functions in the test space @type integer
      %   nTestIC: Number of basis functions in the test space for the initial
      %     condition @type integer

      % save the given values
      obj.pd      = pd;
      obj.nTrial  = nTrial;
      obj.nTest   = nTest;
      obj.nTestIC = nTestIC;
    end
  end

  methods(Abstract)
    % Assemble the spatial mass matrix, that means we evaluate the integral
    % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`.
    %
    % Return values:
    %   Mx: spatial mass matrix @type matrix
    Mx = massMatrix(obj);

    % Assemble the spatial stiffness matrix, that means we evaluate the integral
    % `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j` and `l`.
    %
    % Return values:
    %   Ax: spatial stiffness matrix @type matrix
    Ax = stiffnessMatrix(obj);

    % Evaluate the functional given by a basis function.
    %
    % Evaluates the L2 inner product of a fixed basis function for the given
    % function fun.
    %
    % Parameters:
    %   fun: the argument of the functional @type function_handle
    %   index: index of the basis function to use @type integer
    %
    % Return values:
    %   val: value of the evaluated functional @type double
    val = evaluateFunctional(obj, fun, index);
  end % abstract methods
end % classdef
