% Numerische Lösung mittels Linienmethode / Finite Differenzen für den
% eindimensionalen Fall mit reellwertigem Parameter \omega.
function [output] = onedim(N_Sigmas, plot_all)
if nargin == 0
    N_Sigmas = 3;
    plot_all = 0;
elseif nargin == 1
    plot_all = 0;
end

% Zeit und Ortsintervalle
t_span  = [0, 1];
x_span  = [0, 1];

% Ortsdiskretisierung
x_steps = 40;
x_grid  = linspace(x_span(1), x_span(2), x_steps);

X       = x_grid(2:end-1);
x_h     = x_grid(2) - x_grid(1);
x_N     = x_steps - 2;

% Konstanten
sigma = 1;

% Ansatzfunktionen für Reihenentwicklung
    function phi = ansatz(k)
       if k == 0
           phi = @(x) 0 * x + 1;
       else
           phi = @(x) sqrt(2 / (1 + k^2 * pi^2)) * sin(pi * k * x);
       end
    end

% Anfangsbedingung
    function u = u_0(x)
        u = x.*sin(pi*x);
%         u = exp(-5 * (x - 0.5)^2) - exp(-5/4);
%         u = exp(-5 * (x - 0.5).^2) - exp(-5/4) - 0.5 * (exp(-30 * (x - 0.5).^2) - exp(-30/4));
    end

% Reaktionsterm in Reihe entwicklen
factor = 20;
Phis = {};
for j = 1:N_Sigmas
   Phis{j} = ansatz(j - 1); 
end
Sigmas = factor * (rand(N_Sigmas, 1) - 0.5);
w = generate_w(N_Sigmas, Sigmas);
% plot(x_grid, w(x_grid));
% return;

% Quellterm
    function y = g(t, x)
        y = 0;
    end

% "exakte" Lösung
    function fun = generateSolution(u_0_fun, N)
        %         w = @(x) 0;
        fun = @(t, x) 0;
        for i = 1:N
            coeff = integral(@(x) 2 * u_0_fun(x) .* sin(i * pi *x), 0, 1);
            fun   = @(t, x) fun(t, x) + coeff * sin(i * pi * x) .* exp(-(i^2 * pi^2 * sigma + w(x)) .* t);
        end
    end

% Erzeuge Matrix für Laplace-Operator
tmp = ones(x_N, 1);
M   = (1 / x_h)^2 * spdiags([tmp, -2*tmp, tmp], -1:1, x_N, x_N);

% Erzeuge Matrix für Reaktionsterm
W = zeros(x_N, x_N);
for j = 1:x_N
    for k = 1:N_Sigmas
        W(j, j) = W(j, j) + Sigmas(k) * Phis{k}(j * x_h);
    end
end

% ODE-System
    function dy = ode(t, y)
        Mat = zeros((1 + N_Sigmas) * (x_N));
        Mat(1:x_N, 1:x_N) = sigma * M - W;
        
        for j = 1:N_Sigmas
            rows = (j * x_N +1):((j+1) * x_N);
            colsA = 1:x_N;
            colsB = (j * x_N +1):((j+1) * x_N);
            for k = 1:x_N
                Mat(rows(k), colsA(k)) = - Phis{j}(k * x_h);
            end
            Mat(rows, colsB) = sigma * M - W;
        end
        
        dy  = Mat * y;
        dy(colsA) = dy(colsA) + arrayfun(@(x) g(t, x), x_grid(2:end-1))';
    end

% Anfangswerte
y_0 = zeros((1 + N_Sigmas) * x_N, 1);
y_0(1:x_N) = arrayfun(@(x) u_0(x), x_grid(2:end-1));

% ODE-Solver loslassen
[T, Y] = ode15s(@(t, y) ode(t, y), t_span, y_0);

% Lösung aufteilen in eigentliche ODE und Sensitivity Analysis Teil
[T_mesh, X_mesh] = meshgrid(T, X);
u   =  Y(:, 1:(x_N));

figure(1);
mesh(T_mesh, X_mesh, u');
title('Lösung u');
xlabel('t');
ylabel('x');

if plot_all
    for j = 1:N_Sigmas
        u_sigma  = Y(:, (j*x_N + 1):((j+1)*x_N));
        figure(j + 1);
        mesh(T_mesh, X_mesh, u_sigma');
        title('Sensitivity u_{\sigma_j}');
        xlabel('t');
        ylabel('x');
    end
end

% "exakte" Lösung generieren
% exU   = generateSolution(@(x) u_0(x), 16);
% exUw  = @(t, x) - t .* exU(t, x);
% exUw2 = @(t, x) t.^2 .* exU(t, x);

% Plots
% figure(1);
% mesh(T_mesh, X_mesh, abs(u' - exU(T_mesh, X_mesh)));
% title('Fehler für u');
% xlabel('t');
% ylabel('x');
%
% figure(2);
% mesh(T_mesh, X_mesh, abs(uw' - exUw(T_mesh, X_mesh)));
% title('Fehler für u_\omega');
% xlabel('t');
% ylabel('x');
%
% figure(3);
% mesh(T_mesh, X_mesh, abs(uw2' - exUw2(T_mesh, X_mesh)));
% title('Fehler für u_{\omega\omega}');
% xlabel('t');
% ylabel('x');
end
