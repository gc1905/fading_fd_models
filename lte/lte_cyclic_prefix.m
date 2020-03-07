%% [N_cp_first, N_cp_other, N_sym_slot] = lte_cyclic_prefix(conf, N_fft = 2048)
%%
%% Returns cyclic prefix lengths and symbol count per clot for specified 
%% configuration and FFT size.
%%
%% Arguments:
%%  conf  - cyclic prefix configuration
%%          'none', 'normal', 'extended',
%%          'N_fft'
%%  N_fft - FFT size
%%
%% Returns:
%%  N_cp_first - cyclic prefix length for 1st symbol is slot
%%  N_cp_other - cyclic prefix length for other symbols in slot
%%  N_sym_slot - number of symbols in slot

% Copyright 2017 Grzegorz Cisek (grzegorzcisek@gmail.com)

function [N_cp_first, N_cp_other, N_sym_slot] = lte_cyclic_prefix(conf, N_fft)
  if (nargin < 2)
    N_fft = 2048;
  end

  conf = lower(conf);

  if (strcmp(conf, 'none'))
    N_sym_slot = 7;
    N_cp_first = 0;
    N_cp_other = 0;
  elseif (strcmp(conf, 'normal'))
    N_sym_slot = 7;
    N_cp_first = (10 / 128) * N_fft;
    N_cp_other = ( 9 / 128) * N_fft;
  elseif (strcmp(conf, 'extended'))
    N_sym_slot = 6;
    N_cp_first = (1 / 4) * N_fft;
    N_cp_other = (1 / 4) * N_fft;
  elseif (strcmp(conf, 'n_fft'))
    N_sym_slot = 7;
    N_cp_first = N_fft;
    N_cp_other = N_fft;
  else
    error('invalid cyclic prefix configuration');
  end
end