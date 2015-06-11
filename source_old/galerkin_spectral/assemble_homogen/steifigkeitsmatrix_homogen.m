function [B, B1, B2, B3, B4, B5] = steifigkeitsmatrix_homogen(data, omega, opt, tspan)
% Kurze beschreibung
%
% Die verschiedenen Anteile an der Steifigkeitsmatrix werden durch verschiedene
% Methoden erzeugt und anschließend zusammenaddiert.

    B1 = assemble_first_part(data, opt);
    B2 = assembly_second_part(data, opt);
    % B3 = teil3(data, opt);
    % B4 = teil4(data, opt);
    B3 = 0;
    B4 = 0;
    B5 = assembly_fifth_part(data, opt);

    B = B1 + B2 + B5;
    % B = B1 + B2 + B4 + B5;
end

% FIXME: refactoren und überprüfen
function [M] = assembly_second_part(data, opt)
% Zweite Teilmatrix: \int_Omega \int_I grad u .* grad v dt dx
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
function [B3] = teil3(data, omega, opt)
    % TODO: not implemented!
    error('Not implemented!');
end

% FIXME: refactoren und überprüfen
function [B4] = teil4(data, opt)
    % Zweite Teilmatrix: \int_Omega \int_I grad u .* grad v dt dx
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
function [M] = assembly_fifth_part(data, opt)
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
            Val(ctr) = val;
            ctr = ctr + 1;
        end
    end

    M = sparse(Idy, Idx, Val, opt.dim.Y, opt.dim.X);
end
