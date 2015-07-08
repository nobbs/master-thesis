classdef TemporalAssemblyLinearConstant < TemporalAssemblyAbstract
  % Assemble the temporal mass and stiffness matrices for a basis given through
  % piecewise linear nodal basis functions in the trial temporal part and
  % indicator functions in the test temporal part.

  properties
    tgrid;
  end

  methods

    function obj = TemporalAssemblyLinearConstant(tgrid)
      % Constructor for this assembly class.
      %
      % Parameters:
      %   tgrid: grid for the temporal interval @type vector

      if nargin > 0
        obj.tgrid = tgrid;
      end
    end

    function Mt = massMatrixTrial(obj)
      % checked
      nK = length(obj.tgrid);
      D = diff(obj.tgrid);
      Mt = spdiags([D 0; 0 D]' * [1 2 0; 0 2 1] / 6, [-1 0 1], nK, nK);
    end

    function Mt = massMatrixBoth(obj, tg, ref)
      % checked

      nK = length(tg);

      if ref == 0
        Mt = spdiags(diff(tg)' * [1/2 1/2], 0:1, nK-1, nK);
        return;
      end

      P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
      Mts = obj.massMatrixBoth((P * tg.').', ref - 1);
      Mt = Mts * P;
    end

    function Mt = massMatrixTest(obj, tg, ref)
      % checked

      nK = length(tg) - 1;

      for idx = 1:ref
        P = sparse(interp1(1:(nK + 1), eye((nK + 1)), 1:(1/2):(nK + 1)));
        tg = (P * tg.').';
        nK = length(tg) - 1;
      end

      Mt = sparse(1:nK, 1:nK, abs(diff(tg)));
    end

    function Ct = halfStiffnessMatrix(obj, tg, ref)
      nK = length(tg);

      if ref == 0
        Ct = diff(speye(nK));
        return;
      end

      P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
      Cts = obj.halfStiffnessMatrix((P * tg.').', ref - 1);
      Ct = Cts * P;
    end

    function Mt = massMatrix(obj, kind, useRefinement, spanidx)
      % Assemble the temporal mass Matrix, that means we evaluate the integral
      % `\int_{I} \theta_k(t) \xi_m(t) \diff t` for `k = 1 \dots K` and `m = 1
      % \dots K'` where the temporal basis functions `\theta_k` are piecewise
      % linear nodal basis functions and `\xi_m` are indicatior functions.
      %
      % The parameter nK corresponds to both `K` and `K'`.
      %
      % Parameters:
      %   nK: number of temporal basis functions `K`, `K'` @type integer
      %   kind: which kind of functions to use @type string
      %   spanidx: specify between which time gridpoints the integrals should
      %     be evaluated @type vector
      %
      % Return values:
      %   Mt: temporal mass matrix @type matrix
      %
      % @todo add refinement

      error('Stop!');

      nK = length(obj.tgrid);

      if nargin == 3
        spanidx = [1, nK];
      elseif nargin == 4
        spanidx = [1, nK];
        useRefinement = false;
      end

      switch kind
        case 'both'
          if useRefinement
            P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
            tg = (P * obj.tgrid.').';
            nK = length(tg);
          else
            tg = obj.tgrid;
            nK = length(tg);
          end

          D = diff(tg);

          D(1:spanidx(1) - 1) = 0;
          D(spanidx(2)+1:nK) = 0;

          Mt = spdiags(D.' * [1/2 1/2], [0 1], nK - 1, nK);

          if useRefinement
            Mt = Mt * P;
          end
        case 'trial'
          D = diff(obj.tgrid);
          Mt = spdiags([D 0; 0 D]' * [1 2 0; 0 2 1] / 6, [-1 0 1], nK, nK);
        case 'test'
          nK = length(tg) - 1;

          if useRefinement
            P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
            tg = (P * obj.tgrid.').';
            nK = length(tg);
          else
            tg = obj.tgrid;
            nK = length(tg);
          end

          Mt = spdiags(abs(diff(tgrid)).', 0, nK, nK);
        otherwise
          error();
      end
    end

    function Ct = halfStiffnessMatrix2(obj, useRefinement)
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

      if nargin == 2
        useRefinement = false;
      end

      nK = length(obj.tgrid);

      if useRefinement
        P = sparse(interp1(1:nK, eye(nK), 1:(1/2):nK));
        tg = (P * obj.tgrid.').';
        nK = length(tg);
      else
        tg = obj.tgrid;
        nK = length(tg);
      end

      Ct = diff(speye(nK));

      if useRefinement
        Ct = Ct * P;
      end
    end

    function At = stiffnessMatrix(obj)
      % checked

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

      nK = length(obj.tgrid);
      iD  = 1 ./ diff(obj.tgrid);
      At = spdiags([iD 0; 0 iD]' * [-1 1 0; 0 1 -1], [-1 0 1], nK, nK);
    end

    function et = forwardInitVector(obj)
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

      nK = length(obj.tgrid);
      et = sparse(1, nK);
      et(1) = 1;
    end

    function et = backwardInitVector(obj)
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

      nK = length(obj.tgrid);
      et = sparse(1, nK);
      et(end) = 1;
    end

  end

end
