classdef SpatialAssemblyAbstract < handle
  % Interface for the assembly of spatial structures.
  %
  % Needed for the construction of the space-time-structures, as these are
  % computed through kronecker-products of temporal and spatial matrices and
  % vectors.
  %
  % See also:
  %   TemporalAssemblyAbstract

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

  end % abstract methods

end % classdef
