clear;

N_fft = 512;
f_d = 5;
ITERATIONS = 1;
mp_profile = 'ht';
b = [0:2:16, [20:8:60], [76:16:180], [180:8:200]];

clear FD_MODEL;
FD_MODEL = struct('type', {}, 'num', {}, 'b', {}, 'isi', {}, 'ord', {}, 'label', {});
FD_MODEL(end+1) = struct('type', 'fd', 'num', 1, 'b', 0, 'isi', 0, 'ord', 0, 'label', '\Phi = 0');
FD_MODEL(end+1) = struct('type', 'fd', 'num', 1, 'b', 0, 'isi', 1, 'ord', 0, 'label', 'reduced \Phi');

MARKS = {'b', 'r--', 'k-.s', 'm:x', 'r-h', 'g--*', 'k-<', 'm--*', 'b'};

SNR_accum = zeros(ITERATIONS,length(b),numel(FD_MODEL));

for pb = 1:length(b)
  for i = 1 : ITERATIONS
  	for fidx = 1 : length(FD_MODEL)
  	  FD_MODEL(fidx).b = b(pb);
  	end
    s = fading_model_compare(N_fft, 'normal', 7*20*20, f_d, FD_MODEL, mp_profile, [1,1], 2);
    for idx = 1 : numel(FD_MODEL)
      SNR_accum(i,pb,idx)= s(idx).SNR;
    end
  end
end

constavg = mean(reshape(SNR_accum(:,:,1),1,[]));
constdiff = SNR_accum(:,:,1) - constavg;

SNR_avg = SNR_accum - repmat(constdiff, [1,1,numel(FD_MODEL)]);


hFig = figure;
hold on;

for idx = 1 : numel(FD_MODEL)
  plot(b, 20 * log10(SNR_avg(:,:,idx)), MARKS{idx}, 'DisplayName', FD_MODEL(idx).label);
end

set(hFig, 'pos',[100 100 640 480]);

xlabel('b');
ylabel('SER [dB]');
legend show;
%title(sprintf('N_{FFT} = 256, b = 16, EVA profile.'));
grid on;
hold off;

set(hFig,'Units','Inches');
pos = get(hFig,'Position');
set(hFig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
set(findall(hFig,'-property','FontSize'),'FontSize',11);