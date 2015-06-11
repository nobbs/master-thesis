epsilon = 1;
c = 0.47;
gamma0 = c * (pi / L)^2;
% Cinfty = einbettung();
zeta = 1.2020569;

% K = gamma0 * pi^(2+epsilon) * (4 * Cinfty * L^(3/2) * zeta + pi^(2 + epsilon))^(-1)

% ohne Konstanten Anteil
K = gamma0 * pi^(2 + epsilon) / (sqrt(6) * L^(3/2) * Cinfty * zeta)
% 
% 
% N = 50;
% grid = linspace(0, L, 100);
% 
% for loop = 1:100
%     coeffs = 2 * K * (rand(N, 1) - 0.5);
%     y = coeffs(1) * ones(100, 1);
%     for k = 2:N
%         j = k - 1;
%         y = y + coeffs(k) * (1 / (pi * j))^(1 + epsilon) * sin((pi * j / L) * grid)';
%     end
% %     figure()
%     plot(grid, y);
%     drawnow;
% end
