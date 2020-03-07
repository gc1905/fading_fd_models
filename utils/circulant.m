%C = circulant(A)
% 
% Generates a circulant matrix C from a vector A.
% A may be a column or row vector.
%
% Arguments:
%  A    - first column or row of circulant matrix
%
% Returns:
%  C    - circulant matrix

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [C] = circulant(A)
  if ~isvector(A)
    error('A must be a vector');
  end

  C = zeros(numel(A));
  n = numel(A);

  if iscolumn(A)
    C(:,1) = A;
    for i = 2 : n
      C(:,i) = circshift(A, i-1);
    end 
  else
    C(1,:) = A;
    for i = 2 : n
      C(i,:) = circshift(A, [0, i-1]);
    end 
  end
end