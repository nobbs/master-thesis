classdef AssembleMassMatrix
	% Assemblierung der Masse-Matrix des Galerkin-Verfahrens.
	%
	% Dazu wird, wie bei Galerkin-Verfahren üblich, mit Ansatz- und Testfunktionen
	% gearbeitet. Anders als im Falle der Finite-Elemente-Methode wird hier nicht
	% mit Funktionen mit lokalem Träger sondern globalen Funktionen gearbeitet. Je
	% nach Wahl der Randbedingungen werden andere Basisfunktionen verwendet.
	%
	% Es wird die im Folgenden die LHS der Variationsformulierung `b(u, v) = f(v)`
	% mit `u \in \mathcal{X}_N` und `v = (v_1, v_2) \in \mathcal{Y}_N` bestimmt,
	% wobei
	% ``\begin{aligned}
	%   b(u, v) = &\int_{I} \skp{u_{t}(t)}{v_1(t)}{L_2(\Omega)} + c
	% 	\skp{\nabla u(t)}{\nabla v_1(t)}{L_2(\Omega)} + \mu
	% 	\skp{u(t)}{v_1(t)}{L_2(\Omega)} dt \\&\qquad+
	% 	\skp{u(0)}{v_2}{L_2(\Omega)} + \int_{I} \skp{\omega(x; t)
	% 	u(t)}{v_1(t)}{L_2(\Omega)} dt.
	% \end{aligned}``
	% Der Term mit dem `\omega(x; t)`-Anteil wird separat behandelt, die
	% restlichen Terme werden von links nach rechts durchnummeriert und in jeweils
	% einer eigenen Methode berechnet.
	%
	% Als Ansatz- und Testfunktionen werden momentan Sinus-Funktionen für den Ort
	% und Legendre-Polynome für die Zeit verwendet. Dies entspricht dem Fall
	% homogener Randbedingungen.
	%
	% Der inhomogene Fall wird über eine Fourier-Basis, das heißt, Sinus- und
	% Kosinus-Funktionen abgedeckt.
	% @todo Not yet implemented.

	properties

	end

	methods (Static, Access = public)

		function [B, B1, B2, B3, B4] = assembleWithoutOmega(data, num, dim)
			% Erzeuge Massematrix ohne Omega-Anteil.
			%
			% Parameters:
			%   data: Daten @type struct
			%   num: Anzahl der verwendeten Basisfunktionen @type struct
			%   dim: Dimensionen der endlichdimensionalen Unterräume @type struct
			%
			% Return values:
			%   B: Gesamte Massematrix ohne Omega-Anteil @type sparsematrix
			%   B1: Erster Summand der Massematrix @type sparsematrix
			%   B2: Zweiter Summand der Massematrix @type sparsematrix
			%   B3: Dritter Summand der Massematrix @type sparsematrix
			%   B4: Vierter Summand der Massematrix @type sparsematrix
			%
			% See also:
			%   assembleFirstPart assembleSecondPart assembleThirdPart assembleFourthPart

			% Die einzelnen Summanden berechnen
			B1 = AssembleMassMatrix.assembleFirstPart(num, dim, data.xspan, data.tspan);
			B2 = AssembleMassMatrix.assembleSecondPart(num, dim, data.xspan, data.tspan, data.c_D);
			B3 = AssembleMassMatrix.assembleThirdPart(num, dim, data.xspan, data.tspan, data.mu);
			B4 = AssembleMassMatrix.assembleFourthPart(num, dim, data.xspan);

			% Alles aufsummieren
			B = B1 + B2 + B3 + B4;
		end

	end

	methods (Static, Access = private)

		function M = assembleFirstPart(num, dim, xspan, tspan)
			% Berechnet ersten Summanden der Variationsformulierung.
			%
			% Es wird der erste Summand der Bilinearform der Raum-Zeit-
			% Variationsformulierung für die Ansatz- und Testfunktionen
			% ausgewertet. Konkret wird also das Doppelintegral  ``\int_{I}
			% \int_{\Omega} u_{t}(t, x) v_{1}(t, x) dx dt`` für die Wahl der
			% Basisfunktionen ausgewertet.
			%
			% Für den Fall homogener Randbedingungen werden als Zeit-Funktionen
			% auf ein passendes Zeitintervall `I = [a, b]` geshiftete Legendre-
			% Polynome verwendet, während für die Orts-Funktionen Sinus-
			% Funktionen der Form `s_j(x) = \sin(\frac{\pi j}{L} x)` verwendet
			% werden, wobei `L` die Breite des Orts-Intervalls `\Omega = [0, L]`
			% sei. Dies erlaubt es, obiges Doppelintegral auf die Form
			% ``\int_{T} L'_k L_m dt \int_{\Omega} s_j s_l dx`` zu bringen. Das
			% Orts-Integral nimmt den Wert `\frac{L}{2} \delta_{j l}` an,
			% während das Zeit-Integral sich mittels partieller Integration und
			% der Tatsache, das die Legendre-Polynome eine Basis des
			% Polynomraums darstellen, zu der im Code verwendeten Bedingung
			% umformen, wobei im Falle, dass das Zeit-Integral ungleich Null
			% ist, dieses stets den Wert `\frac{2}{b - a} \cdot 2` annimmt.
			% @todo Prüfen, ob die 2 da wirklich hingehört!
			%
			% Parameters:
			%   num: Anzahl der verwendeten Basisfunktionen @type struct
			%   dim: Dimensionen der endlichdimensionalen Unterräume @type struct
			%   xspan: Grenzen des Orts-Intervalls @type rowvec
			%   tspan: Grenzen des Zeit-Intervalls @type rowvec
			%
			% Return values:
			%   M: Erster Summand der Massematrix @type sparsematrix

			Idx = ones(dim.X, 1);
			Idy = ones(dim.X, 1);
			Val = zeros(dim.X, 1);
			ctr = 1;

			for jdx = 1:min(num.Xj, num.Yl)
				for kdx = 1:num.Xk
					for mdx = 1:kdx
						if kdx > mdx && mod(kdx + mdx, 2) == 1
							Idx(ctr) = (jdx - 1) * num.Xk + kdx;
							Idy(ctr) = (jdx - 1) * num.Ym + mdx;
							Val(ctr) = (xspan(2) / 2) * 2 / (tspan(2) - tspan(1));
							ctr = ctr + 1;
						end
					end
				end
			end

			M = sparse(Idy, Idx, Val, dim.Y, dim.X);
		end

		function M = assembleSecondPart(num, dim, xspan, tspan, diffCoeff)
			% Berechnet den zweiten Summanden der LHS des Variationsproblems.
			%
			% Der zweite Summand (und erster Term der Bilinearform `a`), das
			% heißt ``\int_{I} \int_{\Omega} c \nabla u(t, x) \nabla v_{1}(t,x)
			% dx dt,`` wird ausgewertet.
			%
			% Im homogenen Fall vereinfacht sich dies erneut und kann aufgrund
			% der Orthogonalität der Basisfunktionen ohne numerische Quadratur
			% berechnet werden.
			%
			%  Parameters:
			%   num: Anzahl der verwendeten Basisfunktionen @type struct
			%   dim: Dimensionen der endlichdimensionalen Unterräume @type struct
			%   xspan: Grenzen des Orts-Intervalls @type rowvec
			%   tspan: Grenzen des Zeit-Intervalls @type rowvec
			%   diffCoeff: Der Faktor `c` im Doppelintegral @type double
			%
			% Return values:
			%   M: Zweiter Summand der Massematrix @type sparsematrix

			Idx = ones(dim.X, 1);
			Idy = ones(dim.X, 1);
			Val = zeros(dim.X, 1);
			ctr = 1;

			for jdx = 1:min(num.Xj, num.Yl)
				for kdx = 1:min(num.Xk, num.Ym)
					kk = kdx - 1;

					val = diffCoeff * (tspan(2) - tspan(1)) / (2 * kk + 1) * (pi * jdx)^2 / (2 * xspan(2));

					Idx(ctr) = (jdx - 1) * num.Xk + kdx;
					Idy(ctr) = (jdx - 1) * num.Ym + kdx;
					Val(ctr) = val;
					ctr = ctr + 1;
				end
			end

			M = sparse(Idy, Idx, Val, dim.Y, dim.X);
		end

		function M = assembleThirdPart(num, dim, xspan, tspan, offsetCoeff)
			% Berechnet den dritten Summanden der LHS des Variationsproblems.
			%
			% Der dritte Summand (und zweiter Term der Bilinearform `a`), das
			% heißt ``\int_{I} \int_{\Omega} \mu  u(t, x) v_{1}(t,x) dx dt,``
			% wird ausgewertet.
			%
			% Im homogenen Fall vereinfacht sich dies wieder und kann aufgrund
			% der Orthogonalität der Basisfunktionen ohne numerische Quadratur
			% berechnet werden.
			%
			%  Parameters:
			%   num: Anzahl der verwendeten Basisfunktionen @type struct
			%   dim: Dimensionen der endlichdimensionalen Unterräume @type struct
			%   xspan: Grenzen des Orts-Intervalls @type rowvec
			%   tspan: Grenzen des Zeit-Intervalls @type rowvec
			%   offsetCoeff: Der Faktor `mu` im Doppelintegral @type double
			%
			% Return values:
			%   M: Dritter Summand der Massematrix @type sparsematrix

			L = xspan(2);

			Idx = ones(dim.X, 1);
			Idy = ones(dim.X, 1);
			Val = zeros(dim.X, 1);
			ctr = 1;

			% Da nur im Falle kdx == mdx und jdx == ldx
			for jdx = 1:min(num.Xj, num.Yl)
				ldx = jdx;
				for kdx = 1:min(num.Xk, num.Ym)
					mdx = kdx;
					kk = kdx - 1;
					% mm = mdx - 1;
					% val = (L / 2) * integral(@(t) legendre_dP_shifted(t, kk, data.tspan) .* legendre_P_shifted(t, mm, data.tspan), data.tspan(1), data.tspan(2));
					b = tspan(2);
					a = tspan(1);
					val = offsetCoeff * (b - a) / (2 * kk + 1) * (L / 2);

					if val ~= 0
						x_pos = (jdx - 1) * num.Xk + kdx;
						y_pos = (ldx - 1) * num.Ym + mdx;

						Idx(ctr) = x_pos;
						Idy(ctr) = y_pos;
						Val(ctr) = val;
						ctr = ctr + 1;
					end
					% end
					%   end
					% end
				end
			end

			M = sparse(Idy, Idx, Val, dim.Y, dim.X);
		end

		function M = assembleFourthPart(num, dim, xspan)
			% Berechnet den vierten Summanden der LHS des Variationsproblems.
			%
			% Der vierte Summand, die Anfangsbedingung, das heißt
			% ``\int_{\Omega} c \nabla u(0, x) \nabla v_{2}(x) dx dt,`` wird
			% ausgewertet.
			%
			% Im homogenen Fall vereinfacht sich dies erneut und kann aufgrund
			% der Orthogonalität der Basisfunktionen ohne numerische Quadratur
			% berechnet werden.
			%
			%  Parameters:
			%   num: Anzahl der verwendeten Basisfunktionen @type struct
			%   dim: Dimensionen der endlichdimensionalen Unterräume @type struct
			%   xspan: Grenzen des Orts-Intervalls @type rowvec
			%
			% Return values:
			%   M: Vierter Summand der Massematrix @type sparsematrix

			Idx = ones(dim.X, 1);
			Idy = ones(dim.X, 1);
			Val = zeros(dim.X, 1);
			ctr = 1;

			for jdx = 1:min(num.Xj, num.Yn)
				for kdx = 1:num.Xk
					kk = kdx - 1;

					val = (-1)^(kk) * xspan(2) / 2;

					Idx(ctr) = (jdx - 1) * num.Xk + kdx;
					Idy(ctr) = num.Yl * num.Ym + jdx;
					Val(ctr) = val;
					ctr = ctr + 1;
				end
			end

			M = sparse(Idy, Idx, Val, dim.Y, dim.X);
		end

		function M = assembleOmegaPart(this)
			% @todo Not yet implemented!
			error('Not implemented!');
		end

	end

end
