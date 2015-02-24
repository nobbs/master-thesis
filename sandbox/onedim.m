% Numerische Lösung mittels Linienmethode / Finite Differenzen für den
% eindimensionalen Fall mit reellwertigem Parameter \omega.
function [output] = onedim(plot_all)
    if nargin == 0
        plot_all = 0;
    end

    % Zeit und Ortsintervalle
    t_span  = [0, 1];
    x_span  = [0, 1];
    
    % Ortsdiskretisierung
    x_steps = 50;
    x_grid  = linspace(x_span(1), x_span(2), x_steps);
    
    X       = x_grid(2:end-1);
    x_h     = x_grid(2) - x_grid(1);
    x_N     = x_steps - 2;

    % Konstanten
    sigma = 1;
    
    % Anfangsbedingung
    function u = u_0(x)
        % u = x*sin(pi*x);
        % u = exp(-5 * (x - 0.5)^2) - exp(-5/4);
        u = exp(-5 * (x - 0.5).^2) - exp(-5/4) - 0.5 * (exp(-30 * (x - 0.5).^2) - exp(-30/4));
    end
    
    % Reaktionsterm
    function y = w(x)
        y = 0;
    end
    
    % Quellterm
    function y = g(t, x)
        y = 0;
    end
    
    % "exakte" Lösung
    function fun = generateSolution(u_0_fun, N)
        fun = @(t, x) 0;
        for i = 1:N
            coeff = integral(@(x) 2 * u_0_fun(x) .* sin(i * pi *x), 0, 1);
            fun   = @(t, x) fun(t, x) + coeff * sin(i * pi * x) .* exp(-(i^2 * pi^2 * sigma + w(x)) .* t);
        end
    end
    
    % Erzeuge Massematrix
    tmp = ones(x_N, 1);
    M   = (1 / x_h)^2 * spdiags([tmp, -2*tmp, tmp], -1:1, x_N, x_N);

    tmp = arrayfun(@(x) w(x), x_grid(2:end-1))';
    W   = spdiags(tmp, [0], x_N, x_N);
    
    % ODE-System
    function dy = ode(t, y)
        dy  = zeros(3 * (x_N), 1);
        Mat = [
            sigma * M - W, zeros(x_N), zeros(x_N);
            - eye(x_N), sigma * M - W,  zeros(x_N);
            zeros(x_N), - 2 * eye(x_N), sigma * M - W;
            ];
        dy  = Mat * y + [arrayfun(@(x) g(t, x), x_grid(2:end-1))'; zeros(x_N, 1); zeros(x_N, 1)];
    end

    % Anfangswerte
    y_0 = [arrayfun(@(x) u_0(x), x_grid(2:end-1))'; zeros(x_N, 1); zeros(x_N, 1)];
    
    % ODE-Solver loslassen
    [T, Y] = ode23s(@(t, y) ode(t, y), t_span, y_0);
    
    % Lösung aufteilen in eigentliche ODE und Sensitivity Analysis Teil
    [T_mesh, X_mesh] = meshgrid(T, X);
    u   = Y(:, 1:(x_N));
    uw  = Y(:, (x_N + 1):(2*x_N));
    uw2 = Y(:, (2*x_N + 1):(3*x_N));

    % "exakte" Lösung generieren
    exU   = generateSolution(@(x) u_0(x), 16);
    exUw  = @(t, x) - t .* exU(t, x);
    exUw2 = @(t, x) t.^2 .* exU(t, x);

    % Plots
    figure(1);
    mesh(T_mesh, X_mesh, abs(u' - exU(T_mesh, X_mesh)));
    title('Fehler für u');
    xlabel('t');
    ylabel('x');
    
    figure(2);
    mesh(T_mesh, X_mesh, abs(uw' - exUw(T_mesh, X_mesh)));
    title('Fehler für u_\omega');
    xlabel('t');
    ylabel('x');

    figure(3);
    mesh(T_mesh, X_mesh, abs(uw2' - exUw2(T_mesh, X_mesh)));
    title('Fehler für u_{\omega\omega}');
    xlabel('t');
    ylabel('x');
end
