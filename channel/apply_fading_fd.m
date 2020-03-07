%[iq_fd] = apply_fading_fd(iq, f_d, N_fft, N_cp, model, mpprofile, channel='sos',
%                                mimo=[1,1], isi_en=0, ord=1)
%
% Applies specified channel profile to frequency domain OFDM complex IQ 
% samples using specified channel profile and fading coefficients generators.
%
% Arguments:
%  iq         - complex frequency domain IQ data matrix of size N_sc x N_sym. 
%               Each colomn is a set of mapped subcarriers for single OFDMA 
%               symbol.
%  f_d        - Doppler frequency [Hz]
%  N_fft      - FFT size
%  cp         - cyclic prefix length
%  model      - frequency domain model of Tapped Delay Line:
%               0 - accurate
%               1 - simplified with invariant impulse response
%               2 - linear ICI approximation
%               3 - polynomial ICI approximation (ord = polynomial order)
%  mpprofile  - channel power delay profile string ('epa', 'eva', 'etu') or 2 x L
%               matrix containing multipath tap delays [s] in column 1 and 
%               multipath tap gains [dB] in column 2
%  channel    - method for generation of channel tap coefficients or channel 
%               matrix
%               'zheng' - Rayleigh fading: Zheng and Xiao Sum-of-Sinusoids method
%               'jtc'   - Rayleigh fading: JTC Fader
%               Matrix of size N_sym x (NumOfTaps) containing time varying 
%               channel coefficients may be passed instead.
%  mimo       - vector with MIMO configuration with elements:
%               1 - TX ant num, 2 - RX ant num, 3 - TX ant corr, 4 - RX ant corr
%  b          - band size of channel frequency response matrix reduction
%  isi_en     - if set to non-zero, consider Intersymbol Interference between  
%               OFDMA frames 
%  ord        - polynomial order, valid only in case of model = 3
%
% Returns:
%  iq_fd     - faded iq frequency domain data

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [iq_fade, cfr, cfr_isi] = apply_fading_fd(iq, f_d, frame_cfg, model, mpprofile, alloc, mimo, channel, b, isi_en, ord)
  N_fft = frame_cfg.N_fft;
  f_s = frame_cfg.F_s;
  N_sc = frame_cfg.N_sc;
  sym_in_slot = frame_cfg.N_slot_symbol;

  [cp_first, cp_other] = cyclic_prefix_len(frame_cfg);

  if nargin < 6; alloc = [1,N_sc]; end
  if nargin < 7; mimo = [1,1]; end
  ant_TX = mimo(1);
  ant_RX = mimo(2);
  if ant_TX > 1 || ant_RX > 1
    cor_TX = mimo(3);
    cor_RX = mimo(4);
  end
  MIMO_channels = ant_TX * ant_RX;

  if nargin < 8; channel = 'zheng'; end
  if nargin < 9; b = N_sc; end
  if nargin < 10; isi_en = 1; end
  if nargin < 11; ord = 2; end

  if ~ischar(mpprofile)
    assert(size(mpprofile, 1) == 2, 'mpprofile must be string or L x 2 matrix');
    tds = mpprofile(1,:);
    tgl  = mpprofile(2,:);
  else
    [tds, tgl] = power_delay_profile(mpprofile, 1 / f_s);
  end

  tgl = 10.0 .^ (tgl / 10.0); 

  assert(mod(numel(iq), N_sc) == 0, 'number of elements in iq must be a multiple of subcarriers for specified N_fft');
  assert(length(tds) == length(tgl), 'lengths of tap_gain and tap_delay must be equal');

  N_sym = size(iq,1)*size(iq,2) / N_sc;
  L     = numel(tds);

  % regrid samples
  iq_grid = reshape(iq, [N_sc, N_sym, ant_TX]);

  if ~ischar(channel)
    h = channel;
  else
    if model == 0
    CH_idx = (cp_first+1):(cp_first+N_fft);
    for k = 1:N_sym-1
      if mod(k, sym_in_slot) == 0
        CH_idx(end+1:end+N_fft) = (CH_idx(end)+cp_first+1):(CH_idx(end)+cp_first+N_fft);
      else
        CH_idx(end+1:end+N_fft) = (CH_idx(end)+cp_other+1):(CH_idx(end)+cp_other+N_fft);
      end
    end
    elseif model == 1
      CH_idx = cp_first + N_fft/2;
      for k = 1:N_sym-1
        if mod(k, sym_in_slot) == 0
          CH_idx(end+1) = CH_idx(end) + cp_first + N_fft;
        else
          CH_idx(end+1) = CH_idx(end) + cp_other + N_fft;
        end
      end
    elseif model == 2 || model == 3
      CH_idx = 0;
      for k = 0:N_sym-1
        if mod(k, sym_in_slot) == 0
          CH_idx(end+1) = CH_idx(end) + cp_first + N_fft;
        else
          CH_idx(end+1) = CH_idx(end) + cp_other + N_fft;
        end
      end 
    else
      error('invalid frequency domain model number (0, 1, 2 or 3 supported)');
    end

    if strcmp(channel, 'zheng')
      for m = 1 : MIMO_channels
        for l = 1 : L
          h(:,l,m) = fading_channel_zheng(f_d, f_s, CH_idx);
        end
      end
    elseif strcmp(channel, 'jtc')
      for m = 1 : MIMO_channels
        for l = 1 : L
          h(:,l,m) = fading_channel_jtc(f_d, f_s, CH_idx);
        end
      end
    elseif strcmp(channel, 'idft')
      for m = 1 : MIMO_channels
        for l = 1 : L
          h(:,l,m) = fading_channel_idft(f_d, f_s, CH_idx);
        end
      end
    elseif strcmp(channel, 'none')
      for m = 1 : MIMO_channels
        h(:,:,m) = ones(numel(CH_idx),L);
      end
    else
      error('no such fading channel generation method');
    end
  end

  if (MIMO_channels > 1)
    R = kronecker_correlation_matrix(ant_TX, ant_RX, [cor_TX, cor_RX]);
    C = chol(R);
    
    h_c = zeros(size(h));
    for l = 1 : size(h,2)
      h_c(:,l,:) = (C' * squeeze(h(:,l,:)).').';
    end

    iq_fade_ch = zeros(size(iq_grid,1), size(iq_grid,2), MIMO_channels);
    for m_rx = 1 : ant_RX
      for m_tx = 1 : ant_TX
        ch_id = ((m_rx-1)*ant_RX)+m_tx;
        iq_fade_ch(:,:,ch_id) = sel_tapped_delay_line_fd(model, iq_grid(:,:,m_tx), tds, tgl, h_c(:,:,ch_id), frame_cfg, alloc, b, isi_en, ord);
      end
    end

    iq_fade = zeros(size(iq_grid,1), size(iq_grid,2), ant_RX);
    for m_rx = 1 : ant_RX
      iq_fade(:,:,m_rx) = sum(iq_fade_ch(:,:,(m_rx-1)*ant_TX+1 : m_rx*ant_TX),3);
    end

    h = h_c;
  else
    [iq_fade,cfr,cfr_isi] = sel_tapped_delay_line_fd(model, iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en, ord);
  end
end

function [iq_fade,cfr,cfr_isi] = sel_tapped_delay_line_fd(num, iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en, ord)
  if num == 0
    [iq_fade,cfr,cfr_isi] = tapped_delay_line_fd_m0(iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en);
  elseif num == 1
    [iq_fade,cfr,cfr_isi] = tapped_delay_line_fd_m1(iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en);
  elseif num == 2
    [iq_fade,cfr,cfr_isi] = tapped_delay_line_fd_m2(iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en);
  elseif num == 3
    [iq_fade,cfr,cfr_isi] = tapped_delay_line_fd_m3(iq_grid, tds, tgl, h, frame_cfg, alloc, b, isi_en, ord);
  else
    error('invalid frequency domain model number (0, 1, 2 or 3 supported)');
  end  
end