%[w] = twiddle(k, N)
%
% Returns a matrix of N-point DFT twiddle factors for rotations in matrix k.
%
% Arguments:
%  k    - matrix of integfer rotations 
%  N    - number of DFT bins
%
% Returns:
%  w    - matrix of twiddle factors

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [w] = twiddle(k, N)
  assert(isscalar(N), 'N must be a scalar');
  
  w = exp(-2i * pi * k / N);
end