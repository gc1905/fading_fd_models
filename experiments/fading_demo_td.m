% number of OFDM carriers (128, 256, 512, 1024, 1536, 2048)
N_fft = 2048;

% cyclic prefix configuration
cp = 'normal';

% number of LTE symbols to generate
N_symbols = 7 * 10;

% doppler frequency
f_d = 300;

% channel multipath profile
mp_profile = 'epa';

% if set to 0, send modulation random data instead of Kronecker delta
IMPULSE_RESPONSE = 0;

% modulation order
mod_ord = 4;

frame_cfg = lte_framing_constants(N_fft, cp);
f_s = frame_cfg.F_s;
N_sc = frame_cfg.N_sc;

% % % % % % % % % % % % % 
% TRANSMITTER

if IMPULSE_RESPONSE
  tx_data_mod = ones(N_sc, N_symbols);
else
  tx_data_mod = randiq(mod_ord, [N_sc, N_symbols]);
end

figure(1);
scattermul(reshape(tx_data_mod, N_sc, N_symbols));
title('TX IQ Scatter');
% OFDM modulation
tx_td_data = ofdma_mod(tx_data_mod, frame_cfg);

% % % % % % % % % % % % % 
% WIRELESS CHANNEL
tic;
rx_td_data = apply_fading_td(tx_td_data, f_d, f_s, mp_profile, 'zheng');
toc ;
%rx_td_data = awgn(channeld1, 90);

% % % % % % % % % % % % % 
% RECEIVER

figure(2);
plot([1:length(tx_td_data)], abs(tx_td_data), 'r', [1:length(rx_td_data)], abs(rx_td_data), 'b');
legend('TX', 'RX');
title('Time Domain');

rx_data_mod = ofdma_demod(rx_td_data, frame_cfg, 0);
figure(3);
scattermul(reshape(rx_data_mod, N_sc, N_symbols));
title('RX IQ Scatter');

% time/frequency plot
figure(4);
rx_data_mod_regrid = reshape(rx_data_mod, N_sc, N_symbols);
mesh(abs(rx_data_mod_regrid));
title('TD Fader - channel response');
xlabel('time [symbol]');
ylabel('frequency [subcarrier]');

EVM = evm(rx_data_mod, reshape(tx_data_mod, [], 1));
disp( sprintf('EVM        = %.2f%%', EVM * 100.0) );
power_loss = rms(rx_td_data) / rms(tx_td_data);
disp( sprintf('Power Loss = %.2fdB', 10 * log10(power_loss)) );