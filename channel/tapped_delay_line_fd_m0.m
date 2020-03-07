%[iq_fade] = tapped_delay_line_fd_m0(iq, N_fft, tap_delay, tap_gain, 
%                                          tap_coeff, alloc=[1,N_sc])
%
% Applies specified channel according tapped delay line FIR filter 
% model to IQ data in frequency domain model 0. Channel is considered 
% to be variant during single OFDM time symbol duration - N_fft channel 
% coeffiencient per single path and OFDM symbol.
%
% Arguments:
%  iq        - complex frequency domain iq data matrix of size N_sc x N_sym
%              (each colomn is a set of subcarriers for single OFDM symbol)
%  tap_delay - multipath tap delay vector of size 1 x L (delay unit is 
%              sample index)
%  tap_gain  - multipath tap linear gain vector of size 1 x L
%  tap_coeff - matrix of time domain channel fading coefficients of size 
%              (N_sym*N_fft) x L
%  frame_cfg - frame constants structurte
%  alloc     - two element vector indicating the first and the last 
%              subcarrier allocated for transmission in localized mode. 
%              Default is [1,N_sc] (all subcarriers allocated).
%  b         - band size of channel frequency response matrix reduction
%
% Returns:
%  iq_fade   - faded iq data

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [iq_fade, cfr, cfr_isi] = tapped_delay_line_fd_m0(iq, tap_delay, tap_gain, tap_coeff, frame_cfg, alloc, b, isi_en)
  N_sc = frame_cfg.N_sc;
  N_fft = frame_cfg.N_fft;
  N_sym = size(iq,2);
  [N_cp_first, N_cp_other] = cyclic_prefix_len(frame_cfg);
  N_sym_slot = frame_cfg.N_slot_symbol;
  
  assert(size(iq,1) == N_sc, 'the fist dimension of iq must be equal to N_sc');

  if nargin < 6; alloc = [1, N_sc]; end
  if nargin < 7; b = N_fft; end
  if nargin < 8; isi_en = 1; end

  iq_fade = zeros(N_sc, N_sym);
  L = length(tap_delay);
  alloc_size = alloc(2) - alloc(1) + 1;
  guards = N_fft - N_sc;

  assert(length(tap_delay) == length(tap_gain), 'lengths of tap_gain and tap_delay must be equal');
  assert(all(size(tap_coeff) == [N_sym * N_fft, L]), 'channel must be a matrix of size [N_sym*N_fft,L]');

  Dc = dftmtx(N_fft)' / sqrt(N_fft);
  sc_idx = guards/2+[alloc(1):alloc(2)];

  if nargout > 1
    cfr = zeros(alloc_size,alloc_size,N_sym);
  end
  if nargout > 2
    cfr_isi = zeros(alloc_size,alloc_size,N_sym);
  end

  % construct channel matrix for each symbol
  for u = 1 : N_sym
    if isi_en
      if mod(u, N_sym_slot) == 1
        N_cp = N_cp_first;
      else
        N_cp = N_cp_other;
      end
      [H, B] = cir_mtx(N_fft, tap_delay, tap_gain, tap_coeff(N_fft*(u-1)+1:N_fft*u, :), N_cp);

      Bf = ifft(B, [], 2) * sqrt(N_fft);
      Bf = fftshift(fft(Bf) / sqrt(N_fft));
      Bf = Bf(guards/2+(alloc(1):alloc(2)), guards/2+(alloc(1):alloc(2)));
      Bf = band_mtx_trun(Bf, b);
    else
      H = cir_mtx(N_fft, tap_delay, tap_gain, tap_coeff(N_fft*(u-1)+1:N_fft*u, :), N_fft);
    end

    G = ifft(H, [], 2) * sqrt(N_fft);
    G = fftshift(fft(G) / sqrt(N_fft));
    G = G(guards/2+(alloc(1):alloc(2)), guards/2+(alloc(1):alloc(2)));
    G = band_mtx_trun(G, b);
    
    if nargout > 1
      cfr(:,:,u) = G;
    end
    if nargout > 2
      cfr_isi(:,:,u) = Bf;
    end

    if isi_en == 0 || u == 1
      iq_fade(alloc(1):alloc(2), u) = G * iq(alloc(1):alloc(2), u);
    else
      iq_fade(alloc(1):alloc(2), u) = G * iq(alloc(1):alloc(2), u) + Bf * (iq(alloc(1):alloc(2), u-1));
    end
  end
end