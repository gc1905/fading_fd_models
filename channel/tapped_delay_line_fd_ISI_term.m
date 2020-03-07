% Copyright 2018 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [iq_fade, Bf] = tapped_delay_line_fd_ISI_term(iq, iq_prev, N_fft, N_cp, tap_delay, tap_gain, tap_coeff, alloc, b)
  N_sc = length(iq);

  if nargin < 8; alloc = [1, N_sc]; end
  if nargin < 9; b = N_sc;
  end

  alloc_size = alloc(2) - alloc(1) + 1;
  guards = N_fft - N_sc;

  ISI_idx = find(tap_delay >= N_cp);

  Bf = zeros(alloc_size, alloc_size);
  
  if isempty(ISI_idx)
    iq_fade = zeros(size(iq));
  else
    ISI_tap_pos   = N_fft - tap_delay(ISI_idx) + N_cp + 1;
    ISI_tap_gains = tap_gain(ISI_idx) .* tap_coeff(ISI_idx);
    
    sc_idx = mod(N_fft/2 + guards/2+[alloc(1):alloc(2)], N_fft);

    Cr  = nuidft(ISI_tap_gains, ISI_tap_pos, sc_idx, N_fft) / N_fft; % O{L_ISI * M}
    bfd = nuidft((tap_delay(ISI_idx) - N_cp) .* ISI_tap_gains, ISI_tap_pos, sc_idx, N_fft) / N_fft; % O{L_ISI * M + L_ISI}
    
    EX = 1 ./ (1 - exp(-j * 2 * pi * [1:N_sc-1] / N_fft));
    Pf = exp(- 1j * 2 * pi * N_cp * (guards/2 - 1 + [alloc(1):alloc(2)].') / N_fft);
    
    % O{(2*b+1) * M}
    for n = alloc(1) : alloc(2) % 1 : N_scc
      for m = max(alloc(1), n - b) : min(alloc(2), n + b)
        if n > m
          Bf(n,m) = EX(n-m) * (Cr(m) - Cr(n));
        elseif n < m
          Bf(n,m) = conj(EX(m-n)) * (Cr(m) - Cr(n));
        else
          Bf(n,m) = bfd(m);
        end
      end
    end

    iq_fade = Bf * (iq_prev(alloc(1):alloc(2)) - Pf .* iq(alloc(1):alloc(2))); % O{M + (M+2*b) * M + (M+2*b)}
  end
end