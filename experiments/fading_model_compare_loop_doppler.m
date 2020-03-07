clear;

N_fft = 512;
dopplers = [5:50:1000];
ITERATIONS = 1;
mp_profile = 'eva';
b = 16;
ISI_EN = 0;

clear FD_MODEL;
FD_MODEL = struct('type', {}, 'num', {}, 'b', {}, 'isi', {}, 'ord', {}, 'label', {});
FD_MODEL(end+1) = struct('type', 'fd', 'num', 1, 'b', b, 'isi', ISI_EN, 'ord', 0, 'label', 'Block');
FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', b, 'isi', ISI_EN, 'ord', 1, 'label', 'Poly ord 1 (linear)');
FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', b, 'isi', ISI_EN, 'ord', 2, 'label', 'Poly ord 2');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', b, 'isi', ISI_EN, 'ord', 3, 'label', 'Poly ord 3');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', b, 'isi', ISI_EN, 'ord', 4, 'label', 'Poly ord 4');

% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b',300, 'isi', 1, 'ord', 2, 'label', 'ISI b = 300');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', 24, 'isi', 1, 'ord', 2, 'label', 'ISI b = 24');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', 24, 'isi', 0, 'ord', 2, 'label', 'no ISI b = 24');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 1, 'ord', 2, 'label', 'ISI b = 16');
% FD_MODEL(end+1) = struct('type', 'fd', 'num', 3, 'b', 16, 'isi', 0, 'ord', 2, 'label', 'no ISI b = 16');

[f_s, N_sc, N_rb] = lte_params_from_fft(N_fft);

MARKS = {'b:<', 'r-->', 'k-.s', 'm:x', 'r-h', 'g--*', 'k-<', 'm--*', 'b'};

SNR_accum = zeros(ITERATIONS,length(dopplers), numel(FD_MODEL));

for ph = 1:length(dopplers)
  for i = 1 : ITERATIONS
    s = fading_model_compare(N_fft, 'normal', 7*3, dopplers(ph), FD_MODEL, mp_profile, [1,1], 4);
    for idx = 1 : numel(FD_MODEL)
      SNR_accum(i,ph,idx)= s(idx).SNR;
    end
  end
end

% compute theoretical upper bound of accuracy for given b
SNR_bound = zeros(length(dopplers), 1);
for n = 1 : length(dopplers)
  pwr = jakes_cir_carrier_pwr(N_fft, dopplers(n), f_s, 1, 1 : N_sc/2);
  SNR_bound(n) = 10 * log10( (pwr(1) + 2 * sum(pwr(2:end))) / (2 * sum(pwr(b+2:end))) );
end

hFig = figure;
hold on;

plot(dopplers / 15e3, SNR_bound, MARKS{end}, 'LineWidth', 2, 'DisplayName', 'Theoretical bound [9, Eq. (16)]');

for idx = 1 : numel(FD_MODEL)
  plot(dopplers / 15e3, 20 * log10(SNR_accum(:,:,idx)), MARKS{idx}, 'DisplayName', FD_MODEL(idx).label);
end

axis([0 max(dopplers / 15e3) 10 60]);
set(hFig, 'pos',[100 100 700 480]);

xlabel('Normalized Doppler Frequency');
ylabel('SER [dB]');
legend show;
%title(sprintf('N_{FFT} = 256, b = 16, EVA profile.'));
grid on;
hold off;

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
set(findall(hFig,'-property','FontSize'),'FontSize',11);