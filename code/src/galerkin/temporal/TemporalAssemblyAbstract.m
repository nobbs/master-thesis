classdef TemporalAssemblyAbstract < handle
  % Interface for the assembly of temporal structures.
  %
  % Needed for the construction of the space-time-structures, as these are
  % computed through Kronecker-products of temporal and spatial matrices and
  % vectors.
  %
  % See also:
  %   SpatialAssemblyAbstract

  properties
    % Reference to the problem data object @type ProblemData
    pd;
    % Number of basis functions in the trial space @type integer
    nTrial;
    % Number of basis functions in the test space @type integer
    nTest;
  end

  methods
    function obj = TemporalAssemblyAbstract(pd)
      % Default constructor
      %
      % Parameters:
      %   pd: Reference to the problem data object @type ProblemData

      obj.pd = pd;
    end
  end

  methods(Abstract)
    % Assemble the temporal mass Matrix, that means we evaluate the integral
    % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
    % `m = 1 \dots K'`.
    %
    % Parameters:
    %   space: switch for which space the matrix should be assembled. possible
    %     values are 'trial', 'test' and 'both'. @type string
    %   spanidx: span of the time subinterval to use given in indexes of
    %     pd.tgrid @type vector
    %
    % Return values:
    %   Mt: temporal mass matrix @type matrix
    Mt = massMatrix(obj, space, spanidx);

    % Assemble the temporal "half stiffness" matrix, that means we evaluate the
    % integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
    % `m = 1 \dots K'`.
    %
    % Return values:
    %   Ct: temporal "half stiffness" matrix @type matrix
    Ct = halfStiffnessMatrix(obj);

    % Assemble the temporal stiffness matrix, that means we evaluate the
    % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
    % `k_1, k_2 = 1 \dots K`.
    %
    % Return values:
    %   At: temporal stiffness matrix @type matrix
    At = stiffnessMatrix(obj);

    % Assemble the temporal row vector responsible for the propagation of the
    % initial condition in the case of the forward propagator, that means we
    % evaluate `\theta_k(0)` for `k = 1 \dots K`.
    %
    % Return values:
    %   et: forward propagation vector @type vector
    et = forwardInitVector(obj);

    % Assemble the temporal row vector responsible for the propagation of the
    % initial condition in the case of the forward propagator, that means we
    % evaluate `\theta_k(L)` for `k = 1 \dots K`.
    %
    % Return values:
    %   et: backward propagation vector @type vector
    et = backwardInitVector(obj);
  end % abstract methods

end % classdef
