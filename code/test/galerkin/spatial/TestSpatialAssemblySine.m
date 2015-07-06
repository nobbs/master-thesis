classdef TestSpatialAssemblySine < matlab.unittest.TestCase

  properties
    object;
  end % properties

  properties(TestParameter)
    nX     = struct('small', 5, 'large', 10);
    nY     = struct('small', 5, 'large', 10);
    nC     = struct('small', 5, 'large', 10);
    xwidth = struct('one', 1, 'pi', pi);
  end

  methods(TestMethodSetup)
    function createAssemblyObject(testCase)
      testCase.object = SpatialAssemblySine;
    end
  end

  methods(TestMethodTeardown)
    function deleteAssemblyObject(testCase)
      delete(testCase.object);
    end
  end

  methods(Test)

    function massMatrix(testCase, nX, nY, xwidth)
      testCase.object.xwidth = xwidth;
      actual = testCase.object.massMatrix(nX, nY);

      expected = zeros(nY, nX);
      for xdx = 1:nX
        for ydx = 1:nY
          expected(ydx, xdx) = integral(@(x) testCase.basisFunc(ydx, x, xwidth) .* ...
            testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
        end
      end

      assert(max(max(abs(actual - expected) < 1e-6)));
    end

    function stiffnessMatrix(testCase, nX, nY, xwidth)
      testCase.object.xwidth = xwidth;
      actual = testCase.object.stiffnessMatrix(nX, nY);

      expected = zeros(nY, nX);
      for xdx = 1:nX
        for ydx = 1:nY
          expected(ydx, xdx) = integral(@(x) testCase.basisDerivativeFunc(ydx, x, xwidth) .* testCase.basisDerivativeFunc(xdx, x, xwidth), 0, xwidth);
        end
      end

      assert(max(max(abs(actual - expected) < 1e-6)));
    end

    function fieldDependentSine(testCase, nX, nY, nC, xwidth)
      testCase.object.xwidth = xwidth;
      actual = testCase.object.fieldDependentSine(nX, nY, nC);

      expected = cell(nC, 1);
      for cdx = 1:nC
        expected{cdx} = zeros(nY, nX);
        for xdx = 1:nX
          for ydx = 1:nY
            expected{cdx}(ydx, xdx) = integral(@(x) testCase.sineFieldFunc(cdx, x, xwidth) .* testCase.basisFunc(ydx, x, xwidth) .*testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
          end
        end
      end

      for cdx = 1:nC
        assert(max(max(abs(actual{cdx} - expected{cdx}) < 1e-6)));
      end
    end

    function fieldDependentFourier(testCase, nX, nY, nC, xwidth)
      testCase.object.xwidth = xwidth;
      actual = testCase.object.fieldDependentFourier(nX, nY, nC);

      expected = cell(nC, 1);
      for cdx = 1:nC
        expected{cdx} = zeros(nY, nX);
        for xdx = 1:nX
          for ydx = 1:nY
            expected{cdx}(ydx, xdx) = integral(@(x) testCase.fourierFieldFunc(cdx, x, xwidth) .* testCase.basisFunc(ydx, x, xwidth) .*testCase.basisFunc(xdx, x, xwidth), 0, xwidth);
          end
        end
      end

      for cdx = 1:nC
        assert(max(max(abs(actual{cdx} - expected{cdx}) < 1e-6)));
      end
    end

  end

  methods
    function val = basisFunc(obj, index, x, width)
      val = sin(pi * index * x / width);
    end

    function val = basisDerivativeFunc(obj, index, x, width)
      val = (pi * index / width) * cos(pi * index * x / width);
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
