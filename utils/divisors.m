%d = divisors(n)
% 
% Return a vector of all non-negative dividors of integer n.
%
% Arguments:
%  n    - integer value
%
% Returns:
%  d    - vector of divisors

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [d] = divisors(n)
  k = 1:ceil(n/2);
  d = k( rem(n,k) == 0 );
  d(end+1) = n;
end