%frame_cfg = lte_framing_constants(N_fft, cyclic_prefix_conf)
%
% Generates the OFDM framing constants structure that consolidates
% all the parameters related to OFDM numerology and physical frame 
% structure.
%
% Arguments:
%  N_fft       - FFT size: one of [128 256 512 1024 1536 2048]
%  cyclic_prefix_conf  - cyclic prefix configuration: 'normal' or 'extended'
%
% Returns:
%  frame_cfg - OFDM framing constants structure with the elements:
%              technology - string, fixed to 'LTE'
%              N_fft - FFT size
%              scs - subcarrier spaing in Hz
%              bandwidth - transmission bandwidth in Hz
%              F_s - sampling frequency in Hz
%              N_sc  - number of subcarriers allocated for transmission
%              N_RB - number of PRB allocated for transmission
%              N_sc_RB - number of OFDM subcarriers per PRB (fixed to 12)
%              N_slot_symbol   - number of symbols in a slot
%              N_subframe_slot - number of slots in a subframe
%              N_frame_slot - number of slots in a frame

% Copyright 2018 Grzegorz Cisek (grzegorzcisek@gmail.com)

function frame_cfg = lte_framing_constants(N_fft, cyclic_prefix_conf)
  assert(ismember(N_fft, [128 256 512 1024 1536 2048]), 'N_fft not supported by LTE standard');

  % subcarrier spacing - constant
  scs = 15e3;

  % sampling frequency
  F_s = scs * N_fft;

  % num of subcarriers in Resource Block (RB)
  N_sc_RB = 12;
  
  % determine number of symols in slot
  if strcmp(cyclic_prefix_conf, 'normal')
    N_slot_symbol = 7;
  elseif strcmp(cyclic_prefix_conf, 'extended')
    N_slot_symbol = 6;
  elseif strcmp(cyclic_prefix_conf, 'none')
    N_slot_symbol = 7;
  end

  % determine number of RB and bandwidth based on FFT size
  if N_fft == 128
    N_RB = 6;
    band = 1.4;
  elseif N_fft == 256
    N_RB = 15;
    band = 3;
  elseif N_fft == 512
    N_RB = 25;
    band = 5;
  elseif N_fft == 1024
    N_RB = 50;
    band = 10;
  elseif N_fft == 1536
    N_RB = 75;
    band = 15;
  elseif N_fft == 2048
    N_RB = 100;
    band = 20;
  end

  % other parameters
  N_subframe_slot = 2;
  N_frame_slot = 20;
  N_sc = N_RB * N_sc_RB;

  frame_cfg = struct('technology', 'LTE', 'N_fft', N_fft, 'scs', scs,... 
    'F_s', F_s, 'cp', cyclic_prefix_conf, 'N_sc', N_sc, ...
    'N_RB', N_RB, 'N_sc_RB', N_sc_RB, 'N_slot_symbol', N_slot_symbol,...
    'N_subframe_slot', N_subframe_slot, 'N_frame_slot', N_frame_slot, 'bandwidth', band*1e6);
end