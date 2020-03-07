% compute non-uniform Inverse Digital Fourier Transform
% applies no scaling

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [y] = nuidft(x, xi, yi, N)
  if isrow(x)
  	x = x.';
  end

  W = zeros(length(yi), length(xi));
  for i = 1:length(yi)
    W(i,:) = exp(j * 2 * pi * (yi(i)-1) * (xi-1) / N);
  end

  y = W * x;
end