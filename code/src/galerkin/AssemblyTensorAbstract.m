classdef AssemblyTensorAbstract < handle
  % Assemble the temporal and spatial matrices and vectors which are needed for
  % the assembly of the space-time-matrices via kronecker products.
  %
  % These methods are needed for the assembly of various matrices of the space-
  % time-problem, e.g. the stiffness matrix and the discrete norm matrices.

  methods(Abstract)

    %% Temporal matrices

    % Assemble the temporal mass Matrix, that means we evaluate the integral
    % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
    % `m = 1 \dots K'`.
    %
    % Return values:
    %   MtF: temporal mass matrix @type matrix
    MtF = assembleTemporalMassMatrix(obj);

    % Assemble the temporal "half stiffness" matrix, that means we evaluate the
    % integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
    % `m = 1 \dots K'`.
    %
    % Return values:
    %   CtF: temporal "half stiffness" matrix @type matrix
    CtF = assembleTemporalHalfStiffnessMatrix(obj);

    % Assemble the temporal stiffness matrix, that means we evaluate the
    % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
    % `k_1, k_2 = 1 \dots K`.
    %
    % Return values:
    %   AtF: temporal stiffness matrix @type matrix
    AtF = assembleTemporalStiffnessMatrix(obj);

    % Assemble the temporal row vector responsible for the propagation of the
    % initial condition in the case of the forward propagator, that means we
    % evaluate `\theta_k(0)` for `k = 1 \dots K`.
    %
    % Return values:
    %   etF: forward propagation vector @type vector
    etF = assembleTemporalInitForwardVector(obj);

    % Assemble the temporal row vector responsible for the propagation of the
    % initial condition in the case of the forward propagator, that means we
    % evaluate `\theta_k(L)` for `k = 1 \dots K`.
    %
    % Return values:
    %   etB: backward propagation vector @type vector
    etB = assembleTemporalInitBackwardVector(obj);


    %% Spatial matrices

    % Assemble the spatial mass matrix, that means we evaluate the integral
    % `\int_{\Omega} \sigma_j(x) \sigma_l(x) \diff x` for some `j` and `l`.
    %
    % Return values:
    %   MxF: spatial mass matrix @type matrix
    MxF = assembleSpatialMassMatrix(obj);

    % Assemble the spatial stiffness matrix, that means we evaluate the integral
    % `\int_{\Omega} \sigma'_j(x) \sigma'_l(x) \diff x` for some `j` and `l`.
    %
    % Return values:
    %   AxF: spatial stiffness matrix @type matrix
    AxF = assembleSpatialStiffnessMatrix(obj);

  end % astract methods

end
