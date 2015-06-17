classdef AssemblyAbstract < handle
  % Interface specification for the assembly classes.
  %
  % The main purpose of this interface is to provide a common subset of methods
  % which are essential for the assembly of the stiffness matrix and the load
  % vector.
  %
  % We want to compute the following variational problem
  % ``b(u, v) = f(v), \qquad u \in \mathcal X,~ v = (v_1, v_2) \in \mathcal{Y}``
  % on suitable finite-dimensional subspaces of `\mathcal X` and
  % `\mathcal Y`, where `b` on the left hand side is a bilinear form given by
  % ``b(u, v) = \int_{I} \skp{u_{t}(t)}{v_1(t)}{L_2(\Omega)} + c
  %   \skp{\nabla u(t)}{\nabla v_1(t)}{L_2(\Omega)} + \mu
  %   \skp{u(t)}{v_1(t)}{L_2(\Omega)} + \skp{\omega(x; t)
  %   u(t)}{v_1(t)}{L_2(\Omega)} \diff t + \skp{u(0)}{v_2}{L_2(\Omega)}``
  % and `f` on the ride hand side is a linear functional of the form
  % ``f(v) = \int_{I} \skp{g(t)}{v_1(t)}{L_2(\Omega)} \diff t
  %   + \skp{u_0}{v_2}{L_2(\Omega)}.``
  %
  % As we're trying to reuse most parts of the stiffness matrix for different
  % fields `\omega`, we separate the assembly of the stiffness matrix into
  % the field-independent and the field-dependent parts.

  properties

  end

  methods(Abstract)
    % Assemble the stiffness matrix without the field-dependent term.
    M = assembleStiffnessMatrixWithoutOmega(obj);

    % Assemble only the field-dependent part of the stiffness matrix.
    M = assembleStiffnessMatrixOnlyOmega(obj);

    % Assemble the load vector
    F = assembleRHS(obj);

    % Construct a function handle of the solution
    %
    % Parameters:
    %   solutionCoeffs: solution vector of the linear system @type colvec
    %
    % Return values:
    %   solfun: function handle of the solution function @type function_handle
    solfun = solutionFunctionFromCoeffs(obj, solutionCoeffs);
  end

end
