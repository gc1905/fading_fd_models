%[fd] = ofdma_demod(td, N_fft, cp='none')
%
% Demodulate OFDM modulated data vector td using specified FFT 
% size and cyclic prefix configuration. Number of subcarriers 
% and guardbands is derived form FFT size. Data vector start
% must be aligned to slot start.
% 
% Arguments:
%  td        - time domain data
%  frame_cfg - framing constants structure
%  cppos     - position of cyclic prefix extraction
%
% Returns:
%  fd    - frequency domain carrier data

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [fd] = ofdma_demod(td, frame_cfg, cppos)
  if (nargin < 3)
    cppos = 0;
  end

  [cp_first, cp_other] = cyclic_prefix_len(frame_cfg, 0);

  assert(cppos >= 0 && cppos <= 1.0, 'cppos must be a real number between 0 and 1');

  guards = frame_cfg.N_fft - frame_cfg.N_sc;

  fd = [];
  idx = 1;
  sym = 0;
  cp_len = cp_first;

  while (idx + frame_cfg.N_fft + cp_len - 1 <= length(td))
    cppos_s = round(cp_len * cppos);
    
    tdrs = td(idx+cp_len : idx+cp_len+frame_cfg.N_fft-cppos_s-1);
    tdrs(end+1:end+cppos_s) = td(idx+cp_len-cppos_s:idx+cp_len-1);

    fds = fftshift(fft(tdrs)) / sqrt(frame_cfg.N_fft);

    fd(end+1:end+frame_cfg.N_sc,1) = fds(guards/2+1:guards/2+frame_cfg.N_sc);

    idx = idx + frame_cfg.N_fft + cp_len;
    sym = sym + 1;

    if (mod(sym, frame_cfg.N_slot_symbol) == 0)
      cp_len = cp_first;
    else
      cp_len = cp_other;
    end
  end

end