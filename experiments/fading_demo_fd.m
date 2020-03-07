% number of OFDM carriers (128, 256, 512, 1024, 1536, 2048)
N_fft = 256;

% cyclic prefix configuration
cp = 'normal';

% number of LTE symbols to generate
N_symbols = 7 * 20 ;

% doppler frequency
f_d = 300;

% channel multipath profile
mp_profile = 'etu';

% if set to 0, send modulation random data instead of Kronecker delta
IMPULSE_RESPONSE = 1;

% modulation order
mod_ord = 4;

frame_cfg = lte_framing_constants(N_fft, cp);
f_s = frame_cfg.F_s;
N_sc = frame_cfg.N_sc;

if IMPULSE_RESPONSE
  tx_data_mod = ones(N_sc, N_symbols);
else
  tx_data_mod = randiq(mod_ord, [N_sc, N_symbols]);
end

tic;
sc_data_fade = apply_fading_fd(tx_data_mod, f_d, frame_cfg, 2, mp_profile);
toc;
td_data = ofdma_mod(sc_data_fade, frame_cfg);

figure(1);
plot([1:length(td_data)], abs(td_data));
title('FD Fader - time domain');

figure(2);
mesh(abs(sc_data_fade));
title('FD Fader - channel response');
xlabel('time [symbol]');
ylabel('frequency [subcarrier]');