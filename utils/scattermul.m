%scattermul(x, S=8)
% 
% Displays a scatter plot of a complex signal, assigning different 
% display color to each row of the matrix.
%
% Arguments:
%  x     - complex data matrix
%  S     - size of the marker

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [] = scattermul(x, S)
  if (nargin < 2)
    S = 8;
  end

  [row, col] = size(x);

  colors = jet(row);

  for i = 1 : col
    scatter(real(x(:,i)), imag(x(:,i)), S, colors(i,:));
    hold on;
  end

  xlabel('In-phase');
  ylabel('Quadrature');

  hold off;

end