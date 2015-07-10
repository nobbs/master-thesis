classdef TemporalAssemblyNodal < TemporalAssemblyAbstract
  % Assemble the temporal mass and stiffness matrices for a basis given through
  % piecewise linear nodal basis functions in the trial temporal part and
  % indicator functions in the test temporal part.
  %
  % Refinement of the temporal grid for the test space is also supported and
  % activated by default.

  properties
    % toggle, whether refinement for the test space should be used @type logical
    useRefinement = true;
  end

  methods
    function obj = TemporalAssemblyNodal(pd, useRefinement)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   pd: reference to the problem data object @type ProblemData
      %   useRefinement: toggle whether refinement in the test space should be
      %     used @type logical

      % default values
      if nargin == 1
        useRefinement = true;
      end

      % call the superclass constructor
      obj@TemporalAssemblyAbstract(pd);
      obj.useRefinement = useRefinement;
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
      %     pd.tgrid @type vector
      %
      % Return values:
      %   Mt: temporal mass matrix @type matrix
      %
      % @todo add partial handling!

      % first thing every time
      obj.resetSizes;

      % get the refinement settings
      reflevel = 0;
      if obj.useRefinement
        reflevel = 1;
      end

      % call the right private method
      switch space
        case 'trial'
          Mt = obj.massMatrixTrial;
        case 'test'
          Mt = obj.massMatrixTest(obj.pd.tgrid, reflevel);
        case 'both'
          Mt = obj.massMatrixBoth(obj.pd.tgrid, reflevel);
        otherwise
          error('You specified a non-existing option!');
      end
    end

    function Ct = halfStiffnessMatrix(obj)
      % Assemble the temporal "half stiffness" matrix, that means we evaluate
      % the integral `\int_{I} \theta'_k(t) \xi_m(t) \diff t` for `k = 1 \dots
      % K` and `m = 1 \dots K'`.
      %
      % Return values:
      %   Ct: temporal "half stiffness" matrix @type matrix

      % first thing every time
      obj.resetSizes;

      % get the refinement settings
      reflevel = 0;
      if obj.useRefinement
        reflevel = 1;
      end

      % and call the private method with these settings
      Ct = obj.halfStiffnessMatrixPrivate(obj.pd.tgrid, reflevel);
    end

    function At = stiffnessMatrix(obj)
      % Assemble the temporal stiffness matrix, that means we evaluate the
      % integral `\int_{I} \theta'_{k_1}(t) \theta'_{k_2}(t) \diff t` for
      % `k_1, k_2 = 1 \dots K`.
      %
      % Return values:
      %   At: temporal stiffness matrix @type matrix

      % first thing every time
      obj.resetSizes;

      % and now assemble
      iD = 1 ./ diff(obj.pd.tgrid);
      At = spdiags([iD 0; 0 iD]' * [-1 1 0; 0 1 -1], [-1 0 1], obj.nTrial, ...
                   obj.nTrial);
    end

    function et = forwardInitVector(obj)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(0)` for `k = 1 \dots K`.
      %
      % Return values:
      %   et: forward propagation vector @type vector

      % first thing every time
      obj.resetSizes;

      et = sparse(1, obj.nTrial);
      et(1) = 1;
    end

    function et = backwardInitVector(obj)
      % Assemble the temporal row vector responsible for the propagation of the
      % initial condition in the case of the forward propagator, that means we
      % evaluate `\theta_k(L)` for `k = 1 \dots K`.
      %
      % Return values:
      %   et: backward propagation vector @type vector

      % first thing every time
      obj.resetSizes;

      et = sparse(1, obj.nTrial);
      et(end) = 1;
    end
  end

  methods(Access = 'private')
    function resetSizes(obj)
      % (Re-)Set the number of trial and test basis functions.

      obj.nTrial = length(obj.pd.tgrid);
      obj.nTest  = (1 + obj.useRefinement) * (obj.nTrial - 1);
    end

    function Mt = massMatrixTrial(obj)
      % Assemble the temporal mass matrix for the trial space.
      %
      % Return values:
      %   Mt: mass matrix @type matrix

      D  = diff(obj.pd.tgrid);
      Mt = spdiags([D 0; 0 D]' * [1 2 0; 0 2 1] / 6, [-1 0 1], obj.nTrial, ...
                   obj.nTrial);
    end

    function Mt = massMatrixBoth(obj, tgrid, reflevel)
      % Assemble the temporal mass matrix for the combination of trial and test
      % space.
      %
      % Parameters:
      %   tgrid: temporal grid @type vector
      %   reflevel: refinement level @type integer
      %
      % Return values:
      %   Mt: mass matrix @type matrix

      nK = length(tgrid);
      if reflevel == 0
        Mt = spdiags(diff(tgrid)' * [1/2 1/2], 0:1, nK-1, nK);
        return;
      end

      P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
      Mts = obj.massMatrixBoth((P * tgrid.').', reflevel - 1);
      Mt = Mts * P;
    end

    function Mt = massMatrixTest(~, tgrid, reflevel)
      % Assemble the temporal mass matrix for the test space.
      %
      % Parameters:
      %   tgrid: temporal grid @type vector
      %   reflevel: refinement level @type integer
      %
      % Return values:
      %   Mt: mass matrix @type matrix

      nK = length(tgrid) - 1;
      for idx = 1:reflevel
        P     = sparse(interp1(1:(nK + 1), eye((nK + 1)), 1:(1/2):(nK + 1)));
        tgrid = (P * tgrid.').';
        nK    = length(tgrid) - 1;
      end

      Mt = sparse(1:nK, 1:nK, abs(diff(tgrid)));
    end

    function Ct = halfStiffnessMatrixPrivate(obj, tgrid, reflevel)
      % Assemble the "half stiffness" matrix.
      %
      % This matrix is assembled by integration of the first derivative in the
      % trial space in combination with the test space.
      %
      % Parameters:
      %   tgrid: temporal grid @type vector
      %   reflevel: refinement level @type integer
      %
      % Return values:
      %   Ct: "half stiffness" matrix @type matrix

      nK = length(tgrid);
      if reflevel == 0
        Ct = diff(speye(nK));
        return;
      end

      P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
      Cts = obj.halfStiffnessMatrixPrivate((P * tgrid.').', reflevel - 1);
      Ct = Cts * P;
    end
  end
end
