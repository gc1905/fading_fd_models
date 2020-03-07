% compute non-uniform Digital Fourier Transform
% applies no scaling

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [y] = nudft(x, xi, yi, N, shift)
  if nargin < 5; shift = 0; end

  if isrow(x)
  	x = x.';
  end

  if shift
    yi = yi - N/2;
  end

  W = zeros(length(yi), length(xi));
  for i = 1:length(yi)
    W(i,:) = exp(-j * 2 * pi * (yi(i)-1) * (xi-1) / N);
  end

  y = W * x;
end