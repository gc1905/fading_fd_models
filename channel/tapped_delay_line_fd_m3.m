%[iq_fade] = tapped_delay_line_fd_m3(iq, N_fft, tap_delay, tap_gain, tap_coeff, cp, alloc=[1,N_sc])
%
% Applies specified channel according tapped delay line FIR filter model to IQ data
% in frequency domain model 3. Variations during single OFDM time symbol duration
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

% Copyright 2018 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [iq_fade, cfr, cfr_isi] = tapped_delay_line_fd_m3(iq, tap_delay, tap_gain, tap_coeff, frame_cfg, alloc, b, isi_en, ORD)
  N_sc = frame_cfg.N_sc;
  N_fft = frame_cfg.N_fft;
  N_sym = size(iq,2);
  [N_cp_first, N_cp_other] = cyclic_prefix_len(frame_cfg);
  N_sym_slot = frame_cfg.N_slot_symbol;

  assert(size(iq,1) == N_sc, 'the fist dimension of iq must be equal to N_sc');

  if nargin < 6; alloc = [1, N_sc]; end
  if nargin < 7; b = N_fft; end
  if nargin < 8; isi_en = 1; end
  if nargin < 9; ORD = 2; end

  if ORD == 0
    b = 0;
  end
  if ~ismember(ORD, [0,1,2,3,4])
    warning('specified interpolation order not supported');
  end

  alloc_size = alloc(2) - alloc(1) + 1;

  iq_fade = zeros(N_sc, N_sym);
  L = length(tap_delay);

  guards = N_fft - N_sc;

  assert(length(tap_delay) == length(tap_gain), 'lengths of tap_gain and tap_delay must be equal');
  assert(all(size(tap_coeff) == [N_sym + 1, L]), 'channel must be a matrix of size [1+N_sym,L]');

  q = linspace(- N_fft*0.5 + 0.5, N_fft*0.5 - 0.5, N_fft)';

  q0f = fft(q) / N_fft;
  Q0  = circulant(q0f);
  Q0  = band_mtx_trun(Q0, b);
  Q0r = Q0(guards/2+[alloc(1):alloc(2)], guards/2+[alloc(1):alloc(2)]);
  if ORD > 1
    q1f = fft(q.^2) / N_fft;
    Q1  = circulant(q1f);
    Q1  = band_mtx_trun(Q1, b);
    Q1r = Q1(guards/2+[alloc(1):alloc(2)], guards/2+[alloc(1):alloc(2)]);
  end
  if ORD > 2
    q2f = fft(q.^3) / N_fft;
    Q2  = circulant(q2f);
    Q2  = band_mtx_trun(Q2, b);
    Q2r = Q2(guards/2+[alloc(1):alloc(2)], guards/2+[alloc(1):alloc(2)]);
  end
  if ORD > 3
    q3f = fft(q.^4) / N_fft;
    Q3  = circulant(q3f);
    Q3  = band_mtx_trun(Q3, b);
    Q3r = Q3(guards/2+[alloc(1):alloc(2)], guards/2+[alloc(1):alloc(2)]);
  end

  twidles = exp(-j * 2 * pi * [-N_fft/2:N_fft/2-1]' / N_fft);
  
  if nargout > 1
    cfr = zeros(alloc_size,alloc_size,N_sym);
  end
  if nargout > 2
    cfr_isi = zeros(alloc_size,alloc_size,N_sym);
  end
  
  % precalculate index matrices
  M0 = zeros(1,1,2);
  M0(:,:,1) = inv([- N_fft - N_cp_first]);
  M0(:,:,2) = inv([- N_fft - N_cp_other]);
  M0_idx = ones(N_sym_slot, 1) * 2;
  M0_idx(1) = 1;
  if ORD > 1
    M1 = zeros(2,2,3);
    M1(:,:,1) = inv([ (-N_fft - N_cp_first)^2, (-N_fft - N_cp_first); (N_fft + N_cp_other)^2, (N_fft + N_cp_other)]);
    M1(:,:,2) = inv([ (-N_fft - N_cp_other)^2, (-N_fft - N_cp_other); (N_fft + N_cp_first)^2, (N_fft + N_cp_first)]);
    M1(:,:,3) = inv([ (-N_fft - N_cp_other)^2, (-N_fft - N_cp_other); (N_fft + N_cp_other)^2, (N_fft + N_cp_other)]);    
    M1_idx = ones(N_sym_slot, 1) * 3;
    M1_idx(1) = 1;
    M1_idx(end) = 2;
  end
  if ORD > 2
    M2 = zeros(3,3,4);
    x1 = - 2 * N_fft - N_cp_first - N_cp_other;
    x2 = - N_fft - N_cp_first;
    x3 =   N_fft + N_cp_other;
    M2(:,:,1) = inv([x1^3, x1^2, x1; x2^3, x2^2, x2; x3^3, x3^2, x3]);
    x1 = - 2 * N_fft - 2 * N_cp_other;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_first;    
    M2(:,:,2) = inv([x1^3, x1^2, x1; x2^3, x2^2, x2; x3^3, x3^2, x3]);
    x1 = - 2 * N_fft - N_cp_other - N_cp_first;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_other;
    M2(:,:,3) = inv([x1^3, x1^2, x1; x2^3, x2^2, x2; x3^3, x3^2, x3]);
    x1 = - 2 * N_fft - 2 * N_cp_other;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_other;
    M2(:,:,4) = inv([x1^3, x1^2, x1; x2^3, x2^2, x2; x3^3, x3^2, x3]);
    M2_idx = ones(N_sym_slot, 1) * 4;
    M2_idx(1) = 1;
    M2_idx(2) = 3;
    M2_idx(end) = 2;
  end
  if ORD > 3
    M3 = zeros(4,4,5);
    x1 = - 2 * N_fft - N_cp_first - N_cp_other;
    x2 = - N_fft - N_cp_first;
    x3 =   N_fft + N_cp_other;
    x4 =   2 * N_fft + 2 * N_cp_other; 
    M3(:,:,1) = inv([x1^4, x1^3, x1^2, x1; x2^4, x2^3, x2^2, x2; x3^4, x3^3, x3^2, x3; x4^4, x4^3, x4^2, x4]);
    x1 = - 2 * N_fft - 2 * N_cp_other;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_first;
    x4 =   2 * N_fft + N_cp_other + N_cp_first;  
    M3(:,:,2) = inv([x1^4, x1^3, x1^2, x1; x2^4, x2^3, x2^2, x2; x3^4, x3^3, x3^2, x3; x4^4, x4^3, x4^2, x4]);
    x1 = - 2 * N_fft - N_cp_other - N_cp_first;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_other;
    x4 =   2 * N_fft + 2 * N_cp_other; 
    M3(:,:,3) = inv([x1^4, x1^3, x1^2, x1; x2^4, x2^3, x2^2, x2; x3^4, x3^3, x3^2, x3; x4^4, x4^3, x4^2, x4]);
    x1 = - 2 * N_fft - 2 * N_cp_other;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_other;
    x4 =   2 * N_fft + N_cp_other + N_cp_first;
    M3(:,:,4) = inv([x1^4, x1^3, x1^2, x1; x2^4, x2^3, x2^2, x2; x3^4, x3^3, x3^2, x3; x4^4, x4^3, x4^2, x4]);
    x1 = - 2 * N_fft - 2 * N_cp_other;
    x2 = - N_fft - N_cp_other;
    x3 =   N_fft + N_cp_other;
    x4 =   2 * N_fft + 2 * N_cp_other;
    M3(:,:,5) = inv([x1^4, x1^3, x1^2, x1; x2^4, x2^3, x2^2, x2; x3^4, x3^3, x3^2, x3; x4^4, x4^3, x4^2, x4]);
    M3_idx = ones(N_sym_slot, 1) * 5;
    M3_idx(1) = 1;
    M3_idx(2) = 3;
    M3_idx(end-1) = 4;
    M3_idx(end) = 2;
  end

  for u = 1 : N_sym
    if mod(u, N_sym_slot) == 1
      N_cp = N_cp_first;
    else
      N_cp = N_cp_other;
    end
    
    % Fourier transform of multipath channel impulse response
    G = diag(nudft(tap_gain .* tap_coeff(u+1,:), tap_delay+1, guards/2+[alloc(1):alloc(2)], N_fft, 1));

    %if mod(u, sym_in_slot) == 1
    %  N_cp = N_cp_first;
    %else
    %  N_cp = N_cp_other;
    %end
    %if mod(u+1, sym_in_slot) == 1
    %  N_cp_next = N_cp_first;
    %else
    %  N_cp_next = N_cp_other;
    %end
    
    % calculate polynomial coefficients
    %if u == N_sym || ORD < 2
    %  tap_coeff_delta(1,:) = (tap_coeff(u+1,:) - tap_coeff(u,:)) / (N_fft + N_cp);
    %else
    %  y1 = tap_coeff(u,:) - tap_coeff(u+1,:);
    %  y2 = tap_coeff(u+2,:) - tap_coeff(u+1,:);
    %  x1 = - N_fft - N_cp;
    %  x2 =   N_fft + N_cp_next;
    %  detA = 1 / (x1 * x2 * (x1 - x2));
    %  tap_coeff_delta(2,:) = detA * (  x2   * y1 - x1   * y2);
    %  tap_coeff_delta(1,:) = detA * (- x2^2 * y1 + x1^2 * y2);
    %end

    sym_idx = mod(u-1,N_sym_slot)+1;
    tap_coeff_delta = zeros(ORD, L);

    if u == N_sym || ORD < 2
      y1 = tap_coeff(u,:) - tap_coeff(u+1,:);
      tap_coeff_delta(1,:) = M0(:,:,M0_idx(sym_idx)) * y1;
    elseif u == 1 || ORD < 3
      for l = 1:L
        y1 = tap_coeff(u  ,l) - tap_coeff(u+1,l);
        y2 = tap_coeff(u+2,l) - tap_coeff(u+1,l);
        tap_coeff_delta(2:-1:1,l) = M1(:,:,M1_idx(sym_idx)) * [y1; y2];
      end
    elseif u == N_sym - 1 || ORD < 4
      for l = 1:L
        y1 = tap_coeff(u-1,l) - tap_coeff(u+1,l);
        y2 = tap_coeff(u  ,l) - tap_coeff(u+1,l);
        y3 = tap_coeff(u+2,l) - tap_coeff(u+1,l);
        tap_coeff_delta(3:-1:1,l) = M2(:,:,M2_idx(sym_idx)) * [y1; y2; y3];
      end
    else
      for l = 1:L
        y1 = tap_coeff(u-1,l) - tap_coeff(u+1,l);
        y2 = tap_coeff(u  ,l) - tap_coeff(u+1,l);
        y3 = tap_coeff(u+2,l) - tap_coeff(u+1,l);
        y4 = tap_coeff(u+3,l) - tap_coeff(u+1,l);
        tap_coeff_delta(4:-1:1,l) = M3(:,:,M3_idx(sym_idx)) * [y1; y2; y3; y4];
      end
    end

    % apply linear basis interpolation
    s0 = zeros(alloc_size, 1);
    for l = 1:L
      s0 = s0 + tap_gain(l) * tap_coeff_delta(1,l) * (twidles(guards/2+[alloc(1):alloc(2)]) .^ tap_delay(l));
    end
    G = G + diag(s0) * Q0r;

    % apply quadratic basis interpolation
    if ORD > 1
      s1 = zeros(alloc_size, 1);
      for l = 1:L
        s1 = s1 + tap_gain(l) * tap_coeff_delta(2,l) * (twidles(guards/2+[alloc(1):alloc(2)]) .^ tap_delay(l));
      end
      G = G + diag(s1) * Q1r;
    end
    % apply cubic basis interpolation
    if ORD > 2
      s2 = zeros(alloc_size, 1);
      for l = 1:L
        s2 = s2 + tap_gain(l) * tap_coeff_delta(3,l) * (twidles(guards/2+[alloc(1):alloc(2)]) .^ tap_delay(l));
      end
      G = G + diag(s2) * Q2r;
    end
    % apply 4-th order basis interpolation
    if ORD > 3
      s3 = zeros(alloc_size, 1);
      for l = 1:L
        s3 = s3 + tap_gain(l) * tap_coeff_delta(4,l) * (twidles(guards/2+[alloc(1):alloc(2)]) .^ tap_delay(l));
      end
      G = G + diag(s3) * Q3r;
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