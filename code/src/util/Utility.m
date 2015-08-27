classdef Utility < handle
  % Collection of utility functions as static methods.
  
  methods(Static)
    
    % Conversion of an function evaluated on a equidistant grid to a series
    % expansion and back.
    
    function coeffs = valuesToSineCosine(values)
      % Calculates the fourier coefficients for the given values in terms of
      % sine and cosine coefficients.
      %
      % Warning:
      %   We ignore the zero frequency component as we are working with zero-
      %   mean-value fields.
      %
      % Parameters:
      %   values: values of the periodic function
      %
      % Return values:
      %   coeffs: coefficients of the fourier series. odd indexes are
      %     cosine, even indexes are sine functions.
      
      N = length(values);
      % compute the fft
      fftvalues = fft(values);
      
      % and now get the sine and then the cosine coefficients
      coeffs          = zeros(N, 1);
      coeffs(1:2:end) = 2 * real(fftvalues(2:N/2+1)) / N;
      coeffs(2:2:end) = - 2 * imag(fftvalues(2:N/2+1)) / N;
    end
    
    function values = sineCosineCoeffsToValues(coeffs, xgrid, xwidth)
      
      if nargin == 1
        xwidth = 2 * pi;
        xgrid  = linspace(0, 2*pi, length(coeffs) + 1);
        xgrid  = xgrid(1:end-1);
      end
      
      N = length(coeffs);
      
      values = zeros(size(xgrid, 1), size(xgrid, 2));
      for idx = 1:N
        if mod(idx, 2) == 1
          % odd: cosine
          values = values + coeffs(idx) * cos(pi * (idx + 1) * xgrid / xwidth);
        else
          % even: sine
          values = values + coeffs(idx) * sin(pi * idx * xgrid / xwidth);
        end
      end
    end
    
    function out = dst(y)
      % DST calculate the discreate sign transform of function y.
      
      % Convecrt 1*N vector to N*1
      if size(y,2) > size(y,1)
        y = transpose(y);
      end
      
      N = length(y);
      a = zeros(N,1);
      for m=1:length(a)
        a(m) = ...
          1 / (N - 1) * ...
          ( ...
          sum(y(1:end-1) .* sin(pi * m * (0:N-2)' / (N - 1))) + ...
          sum(y(2:end) .* sin(pi * m * (1:N-1)' / (N - 1))) ...
          );
      end
      out = a;
    end
    
  end
  
  
end
