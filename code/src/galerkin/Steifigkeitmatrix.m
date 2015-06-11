classdef Steifigkeitmatrix
	% Erzeugt die Steifigkeitmatrix für das Galerkin-Verfahen.
	% Doxygen documentation guidlines example class
	%
	% This is the documentation guideline file for class and function
	% documentation. KerMor uses 'doxygen' and a custom tool 'mtoc' to
	% create the documentation from the source files. Doxygen has specific
	% tags to enable easy documentation formatting/layout within the source
	% files. Doxygen commands always begin with an at-character(\@) OR a
	% backslash(\\).

	properties
		data;
		omega;
		opt;
	end

	methods

		function obj = Steifigkeitmatrix(data, omega, opt)
			obj.data = data;
			obj.omega = omega;
			obj.opt = opt;
		end

		function [B, B1, B2, B3, B4, B5] = steifigkeitsmatrix_homogen(this)
			% Kurze beschreibung
			%
			% Die verschiedenen Anteile an der Steifigkeitsmatrix werden durch verschiedene
			% Methoden erzeugt und anschließend zusammenaddiert.

			B1 = this.assemble_first_part();
			B2 = assembly_second_part();
			% B3 = teil3();
			% B4 = teil4();
			B3 = 0;
			B4 = 0;
			B5 = assembly_fifth_part();

			B = B1 + B2 + B5;
			% B = B1 + B2 + B4 + B5;
		end

		function [M] = assemble_first_part(this)
			% Ganz kurze Beschreibung
			%
			% Schauen wir doch mal ob das auch wirklich funktioniert. Hoffen wir mal..
			% bla bla lba bla bie `\int_{0}^{T} f(x) dx`

			Idx = zeros(this.opt.dim.X, 1);
			Idy = zeros(this.opt.dim.X, 1);
			Val = zeros(this.opt.dim.X, 1);
			ctr = 1;

			for jdx = 1:min(this.opt.num.Xj, this.opt.num.Yl)
				for kdx = 1:this.opt.num.Xk
					for mdx = 1:kdx
						if kdx > mdx && mod(kdx + mdx, 2) == 1
							Idx(ctr) = (jdx - 1) * this.opt.num.Xk + kdx;
							Idy(ctr) = (jdx - 1) * this.opt.num.Ym + mdx;
							Val(ctr) = (this.data.xspan(2) / 2) * 2 / (this.data.tspan(2) - this.data.tspan(1));
							ctr = ctr + 1;
						end
					end
				end
			end

			M = sparse(Idy, Idx, Val, this.opt.dim.Y, this.opt.dim.X);

		end

		% FIXME: refactoren und überprüfen
		function [M] = assembly_second_part(this)
			% Zweite Teilmatrix: `\int_Omega \int_I grad u .* grad v dt dx`
			% bla bla lba bla bie `\int_{0}^{T} f(x) dx`
			Idx = ones(opt.dim.X, 1);
			Idy = ones(opt.dim.X, 1);
			Val = zeros(opt.dim.X, 1);
			ctr = 1;

			for jdx = 1:min(opt.num.Xj, opt.num.Yl)
				for kdx = 1:min(opt.num.Xk, opt.num.Ym)
					kk = kdx - 1;

					b = data.tspan(2);
					a = data.tspan(1);
					val = data.c_D * (data.tspan(2) - data.tspan(1)) / (2 * kk + 1) * (pi * jdx)^2 / (2 * data.xspan(2));

					Idx(ctr) = (jdx - 1) * opt.num.Xk + kdx;
					Idy(ctr) = (jdx - 1) * opt.num.Ym + kdx;
					Val(ctr) = val;
					ctr = ctr + 1;
				end
			end

			M = sparse(Idy, Idx, Val, opt.dim.Y, opt.dim.X);

		end

		% FIXME: implementieren
		function [B3] = teil3(this)
			% TODO: not implemented!
			error('Not implemented!');
		end

		% FIXME: refactoren und überprüfen
		function [B4] = teil4(this)
			% Zweite Teilmatrix: `\int_Omega \int_I grad u .* grad v dt dx`
			L = data.xspan(2);

			Idx = zeros(opt.dim.X, 1);
			Idy = zeros(opt.dim.X, 1);
			Val = zeros(opt.dim.X, 1);
			ctr = 1;

			% Da nur im Falle kdx == mdx und jdx == ldx
			for jdx = 1:min(opt.num.Xj, opt.num.Yl)
				ldx = jdx;
				for kdx = 1:min(opt.num.Xk, opt.num.Ym)
					mdx = kdx;
					% for ldx = 1:opt.num.Yl
					%   for mdx = 1:opt.num.Ym
					% FIXME: Erstmal mittels numerischer Quadratur für die
					% unbekannten Teile
					% if jdx == ldx
					kk = kdx - 1;
					% mm = mdx - 1;
					% val = (L / 2) * integral(@(t) legendre_dP_shifted(t, kk, data.tspan) .* legendre_P_shifted(t, mm, data.tspan), data.tspan(1), data.tspan(2));
					b = data.tspan(2);
					a = data.tspan(1);
					val = data.mu * (b - a) / (2 * kk + 1) * (L / 2);

					if val ~= 0
						x_pos = (jdx - 1) * opt.num.Xk + kdx;
						y_pos = (ldx - 1) * opt.num.Ym + mdx;

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

			B4 = sparse(Idy, Idx, Val, opt.dim.Y, opt.dim.X);
		end

		% FIXME: refactoren und überprüfen
		function [M] = assembly_fifth_part(this)
			Idx = zeros(opt.dim.X, 1);
			Idy = zeros(opt.dim.X, 1);
			Val = zeros(opt.dim.X, 1);
			ctr = 1;

			for jdx = 1:min(opt.num.Xj, opt.num.Yn)
				for kdx = 1:opt.num.Xk
					kk = kdx - 1;

					val = (-1)^(kk) * data.xspan(2) / 2;

					Idx(ctr) = (jdx - 1) * opt.num.Xk + kdx;
					Idy(ctr) = opt.num.Yl * opt.num.Ym + jdx;
					Val(ctr) = val;S
					ctr = ctr + 1;
				end
			end

			M = sparse(Idy, Idx, Val, opt.dim.Y, opt.dim.X);
		end

	end

end
