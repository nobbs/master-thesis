classdef GalerkinSolver < handle
  % Löser basierend auf dem Galerkin-Verfahren.
  %
  % @todo Not yet implemented.

  properties
    % Daten des Gebiets

    % Zeitintervall. @type vector
    tspan = [0, 1];
    % Ortsintervall. @type vector
    xspan = [0, 1];


    % Daten des Variationsproblems

    % Anfangsbedingung.
    initialData;
    % Quellterm.
    sourceData;
    % Erstes Feld.
    fieldA;
    % Zweites Feld.
    fieldB;
    % Zeitpunkt, an dem die Felder gewechselt werden. @type double
    fieldBreakpoint = 0.5;
    % Koeffizient vor Laplace-Operator. @type double
    coeffLaplace = 1;
    % Offset-Koeffizient `\mu`. @type double
    coeffOffset = 0;
  end

  properties (Dependent)
    % Check, ob Offset verwendet wird oder nicht. @type logical
    withOffset;
    % Check, ob Anfangsbedingung gesetzt wurde. @type logical
    isInitialDataSet;
    % Check, ob Quellterm gesetzt wurde. @type logical
    isSourceDataSet;
    % Check, ob Felder gesetzt wurden. @type logical
    areFieldsSet;
  end

  methods

    % Konstruktoren

%     function obj = GalerkinSolver()
%       obj;
%     end


    % Getter für Dependent-Properties
    function val = get.withOffset(obj)
      % Prüfe, ob ein Offset bei den Feldern aktiviert ist.
      %
      % Return values:
      %   val: true genau dann wenn, die Felder um eine additive Konstante
      %     verschoben werden @type logical

      if obj.coeffOffset == 0
        val = false;
      else
        val = true;
      end
    end

    function val = get.isInitialDataSet(obj)
      val = true;
    end

    function val = get.isSourceDataSet(obj)
      val = true;
    end

    function val = get.areFieldsSet(obj)
      val = true;
    end


    % Setup-Methoden

    function val = isSetupDone(obj, listFiles)
      % Check, ob das Setup abgeschlossen ist.
      %
      % Geprüft wird, ob die Anfangsbedingungen, Randbedingungen, ein
      % Quellterm et cetera übergeben wurden und alle benötigten Daten für
      % das Lösen der partiellen Differentialgleichung vorliegen.
      %
      % Parameters:
      %   listFiles: Schalter, ob fehlende Daten aufgelistet werden sollen @type logical @default false
      %
      % Return values:
      %   val: Ob Setup abgeschlossen ist oder nicht. @type logical

      if nargin == 1
        listFiles = false;
      end

      if ~obj.isInitialDataSet
        if listFiles
          warning('Anfangsbedingung nicht gesetzt!');
        end
        val = false;
        return;
      end

      if ~obj.isSourceDataSet
        if listFiles
          warning('Quellterm nicht gesetzt!');
        end
        val = false;
        return;
      end

      if ~obj.areFieldsSet
        if listFiles
          warning('Felder nicht gesetzt!');
        end
        val = false;
        return;
      end

      val = true;
    end


    % Setup der Anfangsbedingung

    function setInitialDataFromFunction(obj, ufun)
      % Setze die Anfangbedingung mittels eines Function-Handles.
      %
      % @todo Not yet implemented!
      %
      % Parameters:
      %   ufun: @type function_handle

      error('Not yet implemented!');
    end


    function setInitialDataPointwise(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setInitialDataFromFourierCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setInitialDataFromSineCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end


    % Setup des Quellterms

    function setSourceDataPointwise(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end


    % Setup der Felder

    function setFieldPointwise(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setFieldFromFourierCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end

    function setFieldFromSineCoeffs(obj)
      % @todo Not yet implemented!
      error('Not yet implemented!');
    end


    % Löser

    function [ufun] = solve(obj)
      % @todo Not yet implemented!

      %| @todo Check einbauen!
      % if ~isSetupDone(true)
      %   error('Setup wurde nicht abgeschlossen!');
      % end

      assembly = AssemblySineLegendre();
      assembly.setNumberOfAnsatzFuncs(10, 10);
      assembly.setNumberOfTestFuncsFromAnsatzFuncs();

      assembly.tspan = obj.tspan;
      assembly.xspan = obj.xspan;

      assembly.coeffLaplace = obj.coeffLaplace;
      assembly.coeffOffset = obj.coeffOffset;

      assembly.initialData = @(x) sin(pi * 1 * x / obj.xspan(2));

      LHS = assembly.assembleStiffnessMatrixWithoutOmega();
      RHS = assembly.assembleRHS();

      solfun = @(t, x) assembly.solutionFunctionFromCoeffs(LHS \ RHS, t, x);
      obj.plotSolution(solfun);
    end



    % Plotten und Darstellen

    function plotSolution(obj, solfun)
      % Lösung darstellen als 3D-Plot.
      %
      % Parameters:
      %   solfun: Function-Handle der Lösungsfunktion @type function_handle

      % Gitter erzeugen
      gridt = linspace(obj.tspan(1), obj.tspan(2));
      gridx = linspace(obj.xspan(1), obj.xspan(2));
      [mesht, meshx] = meshgrid(gridt, gridx);

      solution = solfun(mesht, meshx);

      % Plot erstellen
      figure();
      mesh(mesht, meshx, solution);
      title('Lösung');
      xlabel('t');
      ylabel('x');
      zlabel('u');
    end

  end

  methods(Access = Private)

  end

end
