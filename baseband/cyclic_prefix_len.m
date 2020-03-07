function [N_cp_1st, N_cp_other, sym_period] = cyclic_prefix_len(frame_cfg, slot_num)
  if (nargin < 2); slot_num = 0; end
  if strcmp(frame_cfg.technology, 'LTE')
    [N_cp_1st, N_cp_other, sym_period] = lte_cyclic_prefix(frame_cfg.cp, frame_cfg.N_fft);
  else
    error('technology not supported');
  end
end