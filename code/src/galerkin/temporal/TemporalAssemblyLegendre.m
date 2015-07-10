classdef TemporalAssemblyLegendre < TemporalAssemblyAbstract
  % Assemble the temporal mass and stiffness matrices for a basis given through
  % shifted legendre polynomials.
  %
  % Warning:
  %   This class is a leftover from very early stages and only supports the
  %   problem setting with one field, as the switch at the breakpoints isn't
  %   really compatible with global polynomials.

  methods
    function obj = TemporalAssemblyLegendre(pd, nTrial, nTest)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   pd: reference to the problem data object @type ProblemData
      %   nTrial: number of basis functions to use in the trial space
      %     @type integer
      %   nTest: number of basis functions to use in the test space
      %     @type integer

      % call the superclass constructor
      obj@TemporalAssemblyAbstract(pd);
      obj.nTrial = nTrial;
      obj.nTest  = nTest;
    end

    function Mt = massMatrix(obj, space, spanidx)
      % Assemble the temporal mass Matrix, that means we evaluate the integral
      % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and
      % `m = 1 \dots K'`.
      %
      % Parameters:
      %   space: switch for which space the matrix should be assembled. possible
      %     values are 'trial', 'test' and 'both'. @type string
      %   spanidx: span of the time subinterval to use given in indexes of
      %     pd.tgrid. Not supported! @type vector
      %
      % Return values:
      %   Mt: temporal mass matrix @type matrix

      if nargin == 3
        if spanidx(1) ~= 1 || spanidx(end) ~= length(obj.pd.tgrid) - 1
          error('TemporalAssemblyLegendre does not support more than one field!');
        end
      end

      % set the right sizes
      switch space
        case 'trial'
          nY = obj.nTrial;
          nX = obj.nTrial;
        case 'test'
          nY = obj.nTest;
          nX = obj.nTest;
        case 'both'
          nY = obj.nTest;
          nX = obj.nTrial;
        otherwise
          error('You specified a non-existing option!');
      end

      % and compute the matrix
      tmp = (obj.pd.tspan(2) - obj.pd.tspan(1)) ./ (2 * ((1:min(nX, nY)) - 1) + 1);
      Mt  = spdiags(tmp.', 0, nY, nX);
    end

    function Ct = halfStiffnessMatrix(obj)
      % Assemble the temporal "half stiffness" matrix, that means we evaluate
      % the integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots
      % K` and `m = 1 \dots K'`.
      %
      % Return values:
      %   Ct: temporal "half stiffness" matrix @type matrix

      % set the sizes and compute the matrix
      nX = obj.nTrial;
      nY = obj.nTest;
      Ct = spdiags(2 * ones(nY, ceil(nX / 2)), 1:2:nX, nY, nX);
    end

    function At = stiffnessMatrix(obj)
      % Assemble the temporal stiffness matrix, that means we evaluate the
      % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
      % `k_1, k_2 = 1 \dots K`.
      %
      % Return values:
      %   At: temporal stiffness matrix @type matrix

      % set the size
      nX = obj.nTrial;

      % and now assemble the matrix. maybe this loop could be optimized but
      % since this class is deprecated...
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
            intTemporal = 2 * kdx2 * (kdx2 - 1) / (obj.pd.tspan(2) - obj.pd.tspan(1));
          else
            intTemporal = 2 * kdx1 * (kdx1 - 1) / (obj.pd.tspan(2) - obj.pd.tspan(1));
          end
          % save the evaluated integral
          At(kdx1, kdx2) = intTemporal;
        end
      end
    end

    function et = forwardInitVector(obj)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(0)` for `k = 1 \dots K`.
      %
      % Return values:
      %   et: forward propagation vector @type vector

      % same: size and compute
      nX = obj.nTrial;
      et = (-1).^((1:nX) - 1);
    end

    function et = backwardInitVector(obj)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(L)` for `k = 1 \dots K`.
      %
      % Return values:
      %   et: backward propagation vector @type vector

      % same: size and compute
      nX = obj.nTrial;
      et = ones(1, nX);
    end

    function val = basisFunc(obj, index, t)
      % Evaluate the basis function for a given index and t value.
      %
      % Warning:
      %   The enumeration of index starts at 1, that means you'll get the
      %   Legendre polynomial of degree (index - 1)!
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   t: values in which the function should be evaluated @type matrix
      %
      % Return values:
      %   val: values of the basis function in t @type matrix

      val = legendrePolynomial(t, index - 1, obj.tspan);
    end
  end
end
