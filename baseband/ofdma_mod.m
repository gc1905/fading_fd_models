%[td] = ofdma_mod(fd, frame_cfg)
%
% Apply OFDM modulation to input frequency domain data vector
% fd using specified FFT size. Number of subcarriers and 
% guardbands is derived form FFT size.
% 
% Arguments:
%  fd        - frequency domain modulation data
%  frame_cfg - framing constants structure
%
% Returns:
%  td    - time domain samples

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [td] = ofdma_mod(fd, frame_cfg)
  [cp_first, cp_other] = cyclic_prefix_len(frame_cfg);

  assert(mod(numel(fd), frame_cfg.N_sc) == 0, 'number of elements in fd must be a multiple of subcarriers for specified N_fft');

  if (mod(numel(fd), frame_cfg.N_sc * frame_cfg.N_slot_symbol) ~= 0)
    warning('input data vector length is not a multiple of symbol size');
  end

  N_sym = numel(fd) / frame_cfg.N_sc;
  guards = frame_cfg.N_fft - frame_cfg.N_sc;

  fdr = reshape(fd, numel(fd), 1);

  td = [];

  for i = 0 : N_sym - 1
    tds = ifft(ifftshift([zeros(guards/2,1); fdr(1+i*frame_cfg.N_sc:1+(i+1)*frame_cfg.N_sc-1); zeros(guards/2,1)])) * sqrt(frame_cfg.N_fft);

    if (mod(i, frame_cfg.N_slot_symbol) == 0)
      cp_len = cp_first;
    else
      cp_len = cp_other;
    end

    % add cyclic prefix and insert data
    td(end+1:end+cp_len,1) = tds(end-cp_len+1:end);
    td(end+1:end+frame_cfg.N_fft,1)  = tds;
  end
end