span = 1:100;
times = zeros(length(span), 1);

for iter = 1:length(span)
    tic;
    N = span(iter);
    continuity_bound
    times(iter) = toc;
end

plot(span, times);