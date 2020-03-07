%[H, ISI] = cir_mtx(N_fft, tap_delay, tap_gain, tap_coeff, N_cp=N_fft)
% 
% Constructs a Channel Impulse Response matrix of size N_fft. If the size of
% cyclic prefix is specified, ISI term matrix is also constructed.
%
% Arguments:
%  N_fft     - size of FFT
%  tap_delay - multipath tap delay vector of size 1 x L (delay unit is 
%              sample index)
%  tap_gain  - multipath tap linear gain vector of size 1 x L
%  tap_coeff - matrix of time domain channel fading coefficients of size 
%              [N_fft, L] or [1, L]
%  N_cp      - length of cyclic prefix for ISI calculation
%
% Returns:
%  H         - channel impulse response matrix
%  ISI       - ISI term of channel impulse response

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [H, ISI] = cir_mtx(N_fft, tap_delay, tap_gain, tap_coeff, N_cp)
  if nargin < 5
  	N_cp = N_fft;
  end

  L = length(tap_delay);

  assert(length(tap_delay) == length(tap_gain), 'lengths of tap_gain and tap_delay must be equal');

  if all(size(tap_coeff) == [1, L])
    tap_coeff = repmat(tap_coeff, [N_fft, 1]);
  elseif size(tap_coeff) ~= [N_fft, L]
  	error('tap_coeff must contain either N_fft or 1 columns');
  end

  H = zeros(N_fft);
  ISI = zeros(N_fft);

  for n = 1 : N_fft
    for l = 1 : L
      if (tap_delay(l) - n + 1 >= N_cp)
        ISI(n,N_fft-tap_delay(l)+n) = tap_coeff(n,l) * tap_gain(l);
      else
        H(n,mod(n-tap_delay(l)-1,N_fft)+1) = tap_coeff(n,l) * tap_gain(l);
      end
    end
  end
end