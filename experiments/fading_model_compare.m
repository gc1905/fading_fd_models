
function [RESULT] = fading_model_compare(N_fft, cp, N_slots, f_d, fd_models, mp_profile, mimo, IMPULSE_RESPONSE, mod_ord, fir_en)
  % number of OFDM carriers (128, 256, 512, 1024, 1536, 2048)
  if (nargin < 1)
    N_fft = 2048;
  end

  % cyclic prefix configuration
  if nargin < 2
    cp = 'normal';
  end

  frame_cfg = lte_framing_constants(N_fft, cp);
  f_s = frame_cfg.F_s;
  N_sc = frame_cfg.N_sc;

  % number of LTE symbols to generate
  if nargin < 3
    N_slots = 1;
  end
  if strcmp(cp, 'extended')
    N_symbols = 6 * N_slots;
  else
    N_symbols = 7 * N_slots;
  end
    
  % doppler frequency
  if nargin < 4
    f_d = 0;
  end

  % channel multipath profile
  if nargin < 5
    fd_models = struct('type', {}, 'num', {}, 'b', {}, 'isi', {}, 'ord', {}, 'N_fft', {});
    % fd_models(end+1) = struct('type', 'fd', 'num', 1, 'b', 16, 'isi', 0, 'ord', 1, 'N_fft', 0);
    fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 1, 'ord', 1, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 1, 'b', 16, 'isi', 0, 'ord', 0, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 0, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 1, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 2, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 3, 'N_fft', 0);
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 1, 'N_fft', {});
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 2, 'N_fft', {});
    % fd_models(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 3, 'N_fft', {});
    % fd_models(end+1) = struct('type', 'eslami', 'num', 0, 'b', 0, 'isi', 0, 'ord', 0, 'N_fft', 256);
  end

  % channel multipath profile
  if nargin < 6
    mp_profile = 'ht';
  end

  % mimo config
  if nargin < 7
    mimo = [1,1,0.9,0.9];
  end

  % if set to 0, send modulation random data instead of Kronecker delta
  if nargin < 8
    IMPULSE_RESPONSE = 0;
  end

  % modulation order
  if nargin < 9
    mod_ord = 2;
  end

  if 0
    b_fir = radio_filter(81, frame_cfg).';
  else
    b_fir = [0;1;0];
  end
  fir_ntaps = length(b_fir);

  do_plots = (nargout == 0) && 0;
  do_subplot = 0;
  subplot_rows = numel(fd_models) + 1;
  subplot_iter = 0;

  % data generation - common part

  [tap_delay, tap_gain] = power_delay_profile(mp_profile, 1 / f_s);

  if IMPULSE_RESPONSE == 3
    tx_data_mod = exp(-j * rand(N_sc, N_symbols, mimo(1)) * 2 * pi);
  elseif IMPULSE_RESPONSE == 2
    tx_data_mod = (randn(N_sc, N_symbols, mimo(1)) + j * randn(N_sc, N_symbols, mimo(1))) / sqrt(2);
  elseif IMPULSE_RESPONSE == 1
    tx_data_mod = ones(N_sc, N_symbols, mimo(1));
  else 
    tx_data_mod = randiq(mod_ord, [N_sc, N_symbols, mimo(1)]);
  end

  if do_plots && 0
    if do_subplot
      subplot(subplot_rows,2,1);
    else
      figure;
    end
    mesh(abs(tx_data_mod));
    title('Modulation Data');
    xlabel('symbol');
    ylabel('subcarrier');
    zlabel('magnitude');
  end

  % TD FADER CHAIN

  for m = 1 : mimo(1)
    tx_td_fade_td_data(:,m) = ofdma_mod(tx_data_mod(:,:,m), frame_cfg);
  end

  % generate channel coefficients for time domain model (one coefficient per sample)
  N = size(tx_td_fade_td_data, 1);
  L = numel(tap_delay);
  h_fade = zeros(N+1, L, mimo(1)*mimo(2));
  for m = 1 : mimo(1)*mimo(2)
    for l = 1 : L
      h_fade(:,l,m) = fading_channel_zheng(f_d, f_s, [-N_fft/2 1:N]);
    end
  end
  
  % windowing filter
  b_fir(N) = 0;
  b_fir = circshift(b_fir,[-floor(fir_ntaps/2) 0]);
  for m = 1 : mimo(1)
    %tx_td_fade_td_data(:,m) = cconv(b_fir,tx_td_fade_td_data(:,m),N);
    tx_td_fade_td_data(:,m) = tx_td_fade_td_data(:,m);
  end

  h_td_fade = h_fade(2:N+1,:,:);

  tic;
  rx_td_fade_td_data = apply_fading_td(tx_td_fade_td_data, f_d, f_s, [tap_delay; tap_gain], h_td_fade, mimo);
  toc;

  for m = 1 : mimo(2)
    rx_data_mod_regrid(:,:,m) = reshape(ofdma_demod(rx_td_fade_td_data(:,m), frame_cfg), [N_sc, N_symbols]);
  end

  if do_plots && 1
    if do_subplot
      subplot(2,2,1); %subplot(subplot_rows,2,2);
    else
      figure;
    end
    mesh(abs(rx_data_mod_regrid));
    title('Time Domain TDL output');
    xlabel('symbol'); 
    ylabel('subcarrier');
    zlabel('magnitude');
  end

  % FD FADER CHAIN

  % mean CP length per slot
  [cp_1st_sym, cp_next_sym, cp_N_sym] = lte_cyclic_prefix(cp, N_fft);

  RESULT = struct('SNR', {}, 'var', {}, 'sig_pow', {}, 'fderrvec', {});

  for fd_model = fd_models
    if strcmp(fd_model.type, 'fd')
      clear h_fd_fade;
      % resample channel coefficients for FD fader (number of coefficient per symbol)
      if 0 == fd_model.num
        t_idx = 1;
        h_fd_fade = zeros(N_fft * N_symbols, L, mimo(1)*mimo(2));
        for s_idx = 0 : N_symbols-1
          if mod(s_idx, cp_N_sym) == 0
            t_idx = t_idx + cp_1st_sym;
          else
            t_idx = t_idx + cp_next_sym;
          end
          h_fd_fade(s_idx*N_fft+1:(s_idx+1)*N_fft,:,:) = h_td_fade(t_idx:t_idx+N_fft-1,:,:);
          t_idx = t_idx + N_fft;
        end
      elseif 1 == fd_model.num
        t_idx = cp_1st_sym + N_fft / 2;
        h_fd_fade(1,:,:) = h_td_fade(t_idx,:,:);
        for s_idx = 1 : N_symbols-1
          if mod(s_idx, cp_N_sym) == 0
            t_idx = t_idx + N_fft + cp_1st_sym;
          else
            t_idx = t_idx + N_fft + cp_next_sym;
          end
          h_fd_fade(end+1,:,:) = h_td_fade(t_idx,:,:);
        end
      elseif 2 == fd_model.num || 3 == fd_model.num
        t_idx = cp_1st_sym + N_fft / 2;
        h_fd_fade(1,:,:) = h_fade(1,:,:);
        h_fd_fade(2,:,:) = h_td_fade(t_idx,:,:);
        for s_idx = 1 : N_symbols-1
          if mod(s_idx, cp_N_sym) == 0
            t_idx = t_idx + N_fft + cp_1st_sym;
          else
            t_idx = t_idx + N_fft + cp_next_sym;
          end
          h_fd_fade(end+1,:,:) = h_td_fade(t_idx,:,:);
        end
      end

      h_fd_fade = squeeze(h_fd_fade);

      tic;
      tx_fd_fade_data_mod_fade = apply_fading_fd(tx_data_mod, f_d, frame_cfg, fd_model.num, [tap_delay; tap_gain], [1,N_sc], mimo, h_fd_fade, fd_model.b, fd_model.isi, fd_model.ord);
      toc;
      for m = 1 : mimo(1)
        %rx_fd_fade_data(:,m) = cconv(b_fir,ofdma_mod(tx_fd_fade_data_mod_fade(:,:,m), frame_cfg), N);
        rx_fd_fade_data(:,m) = ofdma_mod(tx_fd_fade_data_mod_fade(:,:,m), frame_cfg);
      end
      for m = 1 : mimo(2)
        tx_fd_fade_data_mod_fade(:,:,m) = reshape(ofdma_demod(rx_fd_fade_data(:,m), frame_cfg), [N_sc, N_symbols]);
      end
    elseif strcmp(fd_model.type, 'eslami')
      rx_fd_fade_data = apply_fading_td_fd(tx_td_fade_td_data, f_d, f_s, [tap_delay; tap_gain], h_td_fade, mimo, fd_model.N_fft);
      for m = 1 : mimo(2)
        tx_fd_fade_data_mod_fade(:,:,m) = reshape(ofdma_demod(rx_fd_fade_data(:,m), frame_cfg), [N_sc, N_symbols]);
      end
    else
      error('model type not defined');
    end

    if do_plots && 1
      if do_subplot
        subplot(subplot_rows,2,2+2*subplot_iter+1);
      else
        figure;
      end
      mesh(abs(tx_fd_fade_data_mod_fade));
      plot_title = sprintf('Model %d output', fd_model.num);
      title(plot_title);
      xlabel('symbol');
      ylabel('subcarrier');
      zlabel('magnitude');
    end

    % COMPARE BOTH MODELS

    % time domain
    if do_plots && 1
      figure;
      vrange = min(numel(rx_td_fade_td_data), numel(rx_fd_fade_data));
      plot([1:length(tx_td_fade_td_data)], abs(tx_td_fade_td_data), 'r', [1:length(rx_td_fade_td_data)], abs(rx_td_fade_td_data), 'b', [1:length(rx_fd_fade_data)], abs(rx_fd_fade_data), 'g', [1:vrange], abs(rx_td_fade_td_data(1:vrange) - rx_fd_fade_data(1:vrange)), 'k');
      legend('no fading', 'TD Fader', 'FD Fader', 'difference');
      plot_title = sprintf('Model %d - comparison', fd_model.num);
      title(plot_title);
    end

    % frequency domain error and RMSE
    fd_diff = rx_data_mod_regrid - tx_fd_fade_data_mod_fade;
    SNR = rms(reshape(rx_data_mod_regrid, [], 1)) / rms(reshape(fd_diff, [], 1));
    SNR_dB = 20.0 * log10(SNR);
    VARIANCE = var(reshape(fd_diff, [], 1));
    if do_plots && 0
      if do_subplot
        subplot(2,2,2+subplot_iter); %subplot(subplot_rows,2,2+2*subplot_iter+2);
      else
        figure;
      end
      mesh(abs(fd_diff));
      plot_title = sprintf('Model %d difference (SNR = %.2fdB)', fd_model.num, SNR_dB);
      title(plot_title);
      xlabel('symbol');
      ylabel('subcarrier');
      zlabel('magnitude');
    end

    SNR_sym = rms(rx_data_mod_regrid) ./ rms(fd_diff);
    SNR_sym_dB = 20.0 * log10(SNR_sym);
    if do_plots && 1
      figure(100);
      hold on
      plot(SNR_sym_dB);
      plot_title = sprintf('Model %d SNR per symbol', fd_model.num);
      title(plot_title);
      xlabel('symbol');
      ylabel('SNR [dB]');
    end    

    disp( sprintf('Frequency response SNR = %.2fdB (%f)', SNR_dB, SNR(end)) );
    disp( sprintf('Frequency response error variance  = %.4f', VARIANCE(end)) );

    subplot_iter = subplot_iter + 1;

    RESULT(end+1) = struct('SNR', SNR, 'var', VARIANCE, 'sig_pow', rms(tx_fd_fade_data_mod_fade(:)), 'fderrvec', fd_diff);
  end

  if do_plots
    hold off;
  end
end