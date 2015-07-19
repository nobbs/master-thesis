classdef TestSpatialAssemblyFourier < matlab.unittest.TestCase

  properties
    object;
  end % properties

  properties(TestParameter)
    nX     = struct('small', 5, 'large', 10);
    nC     = struct('small', 5, 'large', 10);
    xwidth = struct('one', 1, 'pi', pi);
  end

  methods(Test)

    function massMatrix(testCase, nX, xwidth)
      pd = ProblemData;
      pd.xspan = [0 xwidth];
      testCase.object = SpatialAssemblyFourier(pd, nX);

      actual = testCase.object.massMatrix;

      expected = zeros(nX, nX);
      for xdx = 1:nX
        for ydx = 1:nX
          expected(ydx, xdx) = integral(@(x) testCase.basisFunc(ydx, x, xwidth) .* ...
            testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
        end
      end

      assert(max(max(abs(actual - expected))) < 1e-6);
    end

    function stiffnessMatrix(testCase, nX, xwidth)
      pd = ProblemData;
      pd.xspan = [0 xwidth];
      testCase.object = SpatialAssemblyFourier(pd, nX);

      actual = testCase.object.stiffnessMatrix;

      expected = zeros(nX, nX);
      for xdx = 1:nX
        for ydx = 1:nX
          expected(ydx, xdx) = integral(@(x) testCase.basisDerivativeFunc(ydx, x, xwidth) .* testCase.basisDerivativeFunc(xdx, x, xwidth), 0, xwidth);
        end
      end

      assert(max(max(abs(actual - expected))) < 1e-6);
    end

    function fieldDependentSine(testCase, nX, nC, xwidth)
      pd = ProblemData;
      pd.xspan = [0 xwidth];
      pd.nC = nC;
      pd.seriesIdx = randi(1000, nC, 1);
      testCase.object = SpatialAssemblyFourier(pd, nX);

      actual = testCase.object.fieldDependentSine;

      expected = cell(nC, 1);
      for cdx = 1:nC
        bdx = pd.seriesIdx(cdx);
        expected{cdx} = zeros(nX, nX);
        for xdx = 1:nX
          for ydx = 1:nX
            expected{cdx}(ydx, xdx) = integral(@(x) testCase.sineFieldFunc(bdx, x, xwidth) .* testCase.basisFunc(ydx, x, xwidth) .*testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
          end
        end
      end

      for cdx = 1:nC
        assert(max(max(abs(actual{cdx} - expected{cdx}))) < 1e-6);
      end
    end

    function fieldDependentFourier(testCase, nX, nC, xwidth)
      pd = ProblemData;
      pd.xspan = [0 xwidth];
      pd.nC = nC;
      pd.seriesIdx = randi(1000, nC, 1);
      testCase.object = SpatialAssemblyFourier(pd, nX);

      actual = testCase.object.fieldDependentFourier;

      expected = cell(nC, 1);
      for cdx = 1:nC
        bdx = pd.seriesIdx(cdx);
        expected{cdx} = zeros(nX, nX);
        for xdx = 1:nX
          for ydx = 1:nX
            expected{cdx}(ydx, xdx) = integral(@(x) testCase.fourierFieldFunc(bdx, x, xwidth) .* testCase.basisFunc(ydx, x, xwidth) .*testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
          end
        end
      end

      for cdx = 1:nC
        assert(max(max(abs(actual{cdx} - expected{cdx}))) < 1e-6);
      end
    end

  end

  methods
    function val = basisFunc(obj, index, x, width)
      % Evaluate the basis function for a given index and x values.
      %
      % Parameters:
      %   index: index of the basis function @type integer
      %   x: values in which the function should be evaluated @type matrix
      %
      % Return values:
      %   val: values of the basis function in x @type matrix

      if index == 1
        val = ones(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = sin(pi * index * x / width);
      else
        val = cos(pi * (index - 1) * x / width);
      end
    end

    function val = basisDerivativeFunc(obj, index, x, width)
      if index == 1
        val = zeros(size(x, 1), size(x, 2));
      elseif mod(index, 2) == 0
        val = (pi * index / width) * cos(pi * index * x / width);
      else
        val = - (pi * (index - 1) / width) * ...
          sin(pi * (index - 1) * x / width);
      end
    end

    function val = sineFieldFunc(obj, index, x, width)
      val = sin(pi * index * x / width);
    end

    function val = fourierFieldFunc(obj, index, x, width)
      if mod(index, 2) == 0
        val = sin(pi * index * x / width);
      else
        val = cos(pi * (index + 1) * x / width);
      end
    end

  end

end
