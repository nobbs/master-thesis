classdef TemporalAssemblyLegendre < TemporalAssemblyAbstract
  % Assemble the temporal mass and stiffness matrices for a basis given through
  % shifted legendre polynomials.

  properties
    % width of the temporal interval @type double
    twidth;
  end

  methods

    function obj = TemporalAssemblyLegendre(twidth)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   twidth: width of the temporal interval @type double

      if nargin > 0
        obj.twidth = twidth;
      end
    end

    function Mt = massMatrix(obj, nX, nY)
      % Assemble the temporal mass Matrix, that means we evaluate the integral
      % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
      % `m = 1 \dots K'` where the temporal basis functions `\theta_k` and
      % `\xi_m` are given as shifted Legendre polynomials.
      %
      % The parameters nX and nY correspond to `K` respectively `K'` and are
      % mainly here, because trial and test space can have a different number
      % of these basis functions.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %   nY: number of temporal basis functions `K'` @type integer
      %
      % Return values:
      %   Mt: temporal mass matrix @type matrix

      tmp = obj.twidth ./ (2 * ((1:min(nX, nY)) - 1) + 1);
      Mt = spdiags(tmp.', 0, nY, nX);
    end

    function Ct = halfStiffnessMatrix(~, nX, nY)
      % Assemble the temporal "half stiffness" matrix, that means we evaluate the
      % integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
      % `m = 1 \dots K'` where the temporal basis functions `\theta_k` and
      % `\xi_m` are defined as shifted Legendre polynomials.
      %
      % The parameters nX and nY correspond to `K` respectively `K'` and are
      % mainly here, because trial and test space can have a different number
      % of these basis functions.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %   nY: number of temporal basis functions `K'` @type integer
      %
      % Return values:
      %   Ct: temporal "half stiffness" matrix @type matrix

      Ct = spdiags(2 * ones(nY, ceil(nX / 2)), 1:2:nX, nY, nX);
    end

    function At = stiffnessMatrix(obj, nX)
      % Assemble the temporal stiffness matrix, that means we evaluate the
      % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
      % `k_1, k_2 = 1 \dots K`, where the temporal basis functions `\theta_k`
      % are defined as shifted Legendre polynomials.
      %
      % As this method is only useful for the assembly of the norm matrices,
      % there's only one parameter nX which corresponds to `K`.
      %
      % Parameters:
      %   nX: number of temporal basis functions `K` @type integer
      %
      % Return values:
      %   At: temporal stiffness matrix @type matrix
      %
      % @todo optimize!

      At = sparse(nX);
      for kdx1 = 1:nX
        % temporal intregal is zero if kdx1 + kdx2 is odd, so we only iterate
        % over the relevant indexes
        if mod(kdx1, 2) == 0
          startKdx2 = 2;
        else
          startKdx2 = 1;
        end
        for kdx2 = startKdx2:2:nX
          % evaluate the temporal integral;
          if kdx1 >= kdx2
            intTemporal = 2 * kdx2 * (kdx2 - 1) / obj.twidth;
          else
            intTemporal = 2 * kdx1 * (kdx1 - 1) / obj.twidth;
          end

          % save the evaluated integral
          At(kdx1, kdx2) = intTemporal;
        end
      end
    end

    function et = forwardInitVector(~, nX)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(0)` for `k = 1 \dots K` with shifted Legendre
      % polynomials `\theta_k`.
      %
      % Parameters:
      %   nX: corresponds to `K` @type integer
      %
      % Return values:
      %   et: forward propagation vector @type vector

      et = (-1).^((1:nX) - 1);
    end

    function et = backwardInitVector(~, nX)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(L)` for `k = 1 \dots K` with shifted Legendre
      % polynomials `\theta_k`.
      %
      % Parameters:
      %   nX: corresponds to `K` @type integer
      %
      % Return values:
      %   et: backward propagation vector @type vector

      et = ones(1, nX);
    end

  end

end
