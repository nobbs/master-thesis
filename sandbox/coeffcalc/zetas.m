%% Beste Exponenten sind anscheinend  p = q = 2
N  = 25;
gr = linspace(2, 10, N);
zp = zeros(N, 1);
zq = zeros(N, 1);

for k = 1:N
    p = gr(k);
    q = p / (p - 1);
    1/p + 1/q
    
    zp(k) = zeta(1 + (1 + ep) * p / 2)^(1/p);
    zq(k) = zeta(1 + (1 + ep) * q / 2)^(1/q);
end

plot(gr, zp, gr, zq, gr, zp .* zq);

%% Bestes beta?
N  = 100;
ep = 0.01;
gr = linspace(1/2 + 0.0001, 3/2 + ep - 0.0001, N);
z1 = zeros(N, 1);
z2 = zeros(N, 1);

for k = 1:N
    bta = gr(k);
    z1(k) = zeta(2*bta)^(1/2);
    z2(k) = zeta(4 + 2 * ep - 2 * bta)^(1/2);
end

plot(gr, z1, gr, z2, gr, z1 .* z2);
