classdef TemporalAssemblyLinearConstant < TemporalAssemblyAbstract
  % Assemble the temporal mass and stiffness matrices for a basis given through
  % piecewise linear nodal basis functions in the ansatz temporal part and
  % indicator functions in the test temporal part.

  methods(Static)

    function Mt = massMatrix(nK, tgrid)
      % Assemble the temporal mass Matrix, that means we evaluate the integral
      % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and `m = 1
      % \dots K'` where the temporal basis functions `\theta_k` are piecewise
      % linear nodal basis functions and `\xi_m` are indicatior functions.
      %
      % The parameter nK corresponds to both `K` and `K'`.
      %
      % Parameters:
      %   nK: number of temporal basis functions `K`, `K'` @type integer
      %   tgrid: temporal grid @type vector
      %
      % Return values:
      %   Mt: temporal mass matrix @type matrix
      %
      % @todo add refinement

      D = diff(tgrid);
      Mt = spdiags(D.' * [1/2 1/2], [0 1], nK - 1, nK);
    end

    function Ct = halfStiffnessMatrix(nK)
      % Assemble the temporal "half stiffness" matrix, that means we evaluate
      % the integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots
      % K` and `m = 1 \dots K'` where the temporal basis functions `\theta_k`
      % are piecewise linear nodal basis functions and `\xi_m` are indicatior
      % functions.
      %
      % The parameter nK corresponds to both `K` and `K'`.
      %
      % Parameters:
      %   nK: number of temporal basis functions `K`, `K'` @type integer
      %
      % Return values:
      %   Ct: temporal "half stiffness" matrix @type matrix
      %
      % @todo add refinement

      Ct = diff(speye(nK));
    end

    function At = stiffnessMatrix(nK, tgrid)
      % Assemble the temporal stiffness matrix, that means we evaluate the
      % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for `k_1,
      % k_2 = 1 \dots K`, where the temporal basis functions `\theta_k` are
      % piecewise linear nodal basis functions.
      %
      % As this method is only useful for the assembly of the norm matrices,
      % there's only one parameter nX which corresponds to `K`.
      %
      % Parameters:
      %   nK: number of temporal basis functions `K` @type integer
      %
      % Return values:
      %   At: temporal stiffness matrix @type matrix

      iD  = 1 ./ diff(tgrid);
      At = spdiags([iD 0; 0 iD]' * [-1 1 0; 0 1 -1], [-1 0 1], nK, nK);
    end

    function et = forwardInitVector(nK)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator that means we
      % evaluate `\theta_k(0)` for `k = 1 \dots K` with piecewise linear nodal
      % basis functions `\theta_k`.
      %
      % Parameters:
      %   nK: corresponds to `K` @type integer
      %
      % Return values:
      %   et: forward propagation vector @type vector

      et = zeros(1, nK);
      et(1) = 1;
    end

    function et = backwardInitVector(nK)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(L)` for `k = 1 \dots K` with piecewise linear nodal
      % basis functions `\theta_k`.
      %
      % Parameters:
      %   nK: corresponds to `K` @type integer
      %
      % Return values:
      %   et: backward propagation vector @type vector

      et = zeros(1, nK);
      et(end) = 1;
    end

  end

end
