% calculates B = dftmtx(N_fft) * A * dftmtx(N_fft)'

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function B = mtx_dft(A)
  assert(size(A,1) == size(A,2), 'A must be square matrix');
  N_fft = size(A,1);

  B = ifft(A, [], 2) * sqrt(N_fft);
  B = fft(B) / sqrt(N_fft);
end