%[iq_fade] = tapped_delay_line_fd_m2(iq, N_fft, tap_delay, tap_gain, tap_coeff, cp, alloc=[1,N_sc])
%
% Applies specified channel according tapped delay line FIR filter model to IQ data
% in frequency domain model 2. Variations during single OFDM time symbol duration
% are linearly approximated - one channel coeffiencient per single path and OFDM symbol.
%
% Arguments:
%  iq        - complex frequency domain iq data matrix of size N_sc x N_sym (each
%              colomn is a set of subcarriers for single OFDM symbol)
%  tap_delay - multipath tap delay vector of size 1 x L (delay unit is sample index)
%  tap_gain  - multipath tap linear gain vector of size 1 x L
%  tap_coeff - matrix of time domain channel fading coefficients of size (1+N_sym) x L
%  frame_cfg - frame constants structurte
%  alloc     - two element vector indicating the first and the olast subcarrier allocated for
%              transmission in localized mode. Default is [1,N_sc] (all subcarriers allocated).
%  b         - band size of channel frequency response matrix reduction
%
% Returns:
%  iq_fade   - faded iq data

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [iq_fade, cfr, cfr_isi] = tapped_delay_line_fd_m2(iq, tap_delay, tap_gain, tap_coeff, frame_cfg, alloc, b, isi_en)
  N_sc = frame_cfg.N_sc;
  N_fft = frame_cfg.N_fft;
  N_sym = size(iq,2);
  [N_cp_first, N_cp_other] = cyclic_prefix_len(frame_cfg);
  N_sym_slot = frame_cfg.N_slot_symbol;

  assert(size(iq,1) == N_sc, 'the fist dimension of iq must be equal to N_sc');

  if nargin < 6; alloc = [1, N_sc]; end
  if nargin < 7; b = N_fft; end
  if nargin < 8; isi_en = 1; end
  
  alloc_size = alloc(2) - alloc(1) + 1;
  iq_fade = zeros(N_sc, N_sym);
  L = length(tap_delay);
  guards = N_fft - N_sc;

  assert(length(tap_delay) == length(tap_gain), 'lengths of tap_gain and tap_delay must be equal');
  assert(all(size(tap_coeff) == [N_sym + 1, L]), 'channel must be a matrix of size [1+N_sym,L]');

  LIN = circulant(fft(linspace(N_fft*0.5 - 0.5,-N_fft*0.5 + 0.5,N_fft).') / N_fft);
  twidles = exp(-j * 2 * pi * [-N_fft/2:N_fft/2-1]' / N_fft);
  LINr = band_mtx_trun(LIN(guards/2+[alloc(1):alloc(2)],guards/2+[alloc(1):alloc(2)]), b);

  if nargout > 1
    cfr = zeros(alloc_size,alloc_size,N_sym);
  end
  if nargout > 2
    cfr_isi = zeros(alloc_size,alloc_size,N_sym);
  end

  for u = 1 : N_sym
    G = diag(nudft(tap_gain .* tap_coeff(u+1,:), tap_delay+1, guards/2+[alloc(1):alloc(2)], N_fft, 1));

    if mod(u, N_sym_slot) == 1
      N_cp = N_cp_first;
    else
      N_cp = N_cp_other;
    end
    tap_coeff_delta = (tap_coeff(u,:) - tap_coeff(u+1,:)) / (N_fft + N_cp);

    for l = 1:L
      G = G + bsxfun(@times, LINr, tap_gain(l) * tap_coeff_delta(l) * (twidles(guards/2+[alloc(1):alloc(2)]) .^ tap_delay(l)));
    end

    if nargout > 1
      cfr(:,:,u) = G;
    end

    % frequency domain convolution
    iq_fade(alloc(1):alloc(2),u) = G * iq(alloc(1):alloc(2),u);

    % ISI effect
    if isi_en
      if u ~= 1
        [iq_isi, Bf] = tapped_delay_line_fd_ISI_term(iq(:,u), iq(:,u-1), N_fft, N_cp, tap_delay, tap_gain, tap_coeff(u,:), alloc, b);
      else
        [iq_isi, Bf] = tapped_delay_line_fd_ISI_term(iq(:,u), zeros(N_sc,1), N_fft, N_cp, tap_delay, tap_gain, tap_coeff(u,:), alloc, b);
      end
      iq_fade(alloc(1):alloc(2),u) = iq_fade(alloc(1):alloc(2),u) + iq_isi;
      if nargout > 2
        cfr_isi(:,:,u) = Bf;
      end
    end
  end
end