function [w] = generate_w( N, sigma )
    w = @(x) 0*x + sigma(1);
    
    for j = 1:(N-1)
       w = @(x) w(x) + sqrt(2 / (1+ j^2 * pi^2)) * sigma(j + 1) * sin(pi * j * x); 
    end
end

