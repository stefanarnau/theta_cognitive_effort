clear all;

% sampling rate
srate = 200;

% Define frequencies
freqs = 1 : 0.2 : 30;

% Number of trials
n_epochs = 200;

% time vector
time_vec = -0.8 : (1 / srate) : 1.8;

% Create simulated data (all amplitude 1)
d = zeros(length(freqs), length(time_vec), n_epochs / 2);
for e = 1 : n_epochs / 2
    for f = 1 : length(freqs)
        phase_offset = -pi + rand(1) * (pi - (-pi));
        d(f, :, e) = sin(2 * pi * freqs(f) * time_vec + phase_offset);
    end
end

% Copy data for two conditions
d = cat(3, d, d);

% Create tf-datasets with baseline differences
idx_1 = 1 : n_epochs / 2;
idx_2 = (n_epochs / 2) + 1 : n_epochs;
tf_data = d;
tf_data(:, time_vec < 0, idx_1) = tf_data(:, time_vec < 0, idx_1) .* 0.5;
tf_data(:, time_vec < 0, idx_2) = tf_data(:, time_vec < 0, idx_2) .* 1.5;

% Sum frequencies
tf_data = squeeze(sum(tf_data, 1));

% Set complex Morlet wavelet parameters
n_frq = 20;
frqrange = [2, 20];
tfres_range = [600, 300];

% Set wavelet time
wtime = -2 : 1 / srate : 2;

% Create wavelet frequencies and tapering Gaussian widths in temporal domain
tf_freqs = logspace(log10(frqrange(1)), log10(frqrange(2)), n_frq);
fwhmTs = logspace(log10(tfres_range(1)), log10(tfres_range(2)), n_frq);

% Init matrices for wavelets
cmw = zeros(length(tf_freqs), length(wtime));
cmwX = zeros(length(tf_freqs), length(wtime));

% Create the wavelets
for frq = 1 : length(tf_freqs)

    % Create wavelet with tapering gaussian corresponding to desired width in temporal domain
    cmw(frq, :) = exp(2 * 1i * pi * tf_freqs(frq) .* wtime) .* exp((-4 * log(2) * wtime.^2) ./ (fwhmTs(frq) / 1000)^2);

    % Normalize wavelet
    cmw(frq, :) = cmw(frq, :) ./ max(cmw(frq, :));

end

% Define time window of analysis
prune_times = [-0.5, 1]; 
tf_times = time_vec(dsearchn(time_vec', prune_times(1)) : dsearchn(time_vec', prune_times(2)));

% Init tf matrix
powcube = NaN(length(tf_freqs), length(time_vec), n_epochs);

% convolution length
convlen = length(time_vec) * n_epochs + size(cmw, 2) - 1;

% cmw to freq domain and scale
cmwX = zeros(length(tf_freqs), convlen);
for f = 1 : length(tf_freqs)
    cmwX(f, :) = fft(cmw(f, :), convlen);
    cmwX(f, :) = cmwX(f, :) ./ max(cmwX(f, :));
end

% Get TF-power
tmp = fft(reshape(tf_data, 1, []), convlen);
for f = 1 : length(tf_freqs)
    as = ifft(cmwX(f, :) .* tmp); 
    as = as(((size(cmw, 2) - 1) / 2) + 1 : end - ((size(cmw, 2) - 1) / 2));
    as = reshape(as, length(time_vec), n_epochs);
    powcube(f, :, :) = abs(as) .^ 2;   
end
   
% Cut edges
powcube = powcube(:, dsearchn(time_vec', prune_times(1)) : dsearchn(time_vec', prune_times(2)), :);

% Get condition general baseline values
ersp_bl = [-0.5, -0.2];
tmp = squeeze(mean(powcube, 3));
[~, blidx1] = min(abs(tf_times - ersp_bl(1)));
[~, blidx2] = min(abs(tf_times - ersp_bl(2)));
blvals = squeeze(mean(tmp(:, blidx1 : blidx2), 2));

% Calculate ersp
raw_1 = squeeze(mean(powcube(:, :, idx_1), 3));
raw_2 = squeeze(mean(powcube(:, :, idx_2), 3));
ersp_1 = 10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, idx_1), 3)), blvals));
ersp_2 = 10 * log10(bsxfun(@rdivide, squeeze(mean(powcube(:, :, idx_2), 3)), blvals));

% Plot ersp
figure()

subplot(2, 3, 1)
contourf(tf_times, tf_freqs, raw_1, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('power - condition 1')

subplot(2, 3, 2)
contourf(tf_times, tf_freqs, raw_2, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('power - condition 2')

subplot(2, 3, 3)
contourf(tf_times, tf_freqs, raw_1 - raw_2, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('power - condition difference')

subplot(2, 3, 4)
contourf(tf_times, tf_freqs, ersp_1, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('normalized - condition 1')

subplot(2, 3, 5)
contourf(tf_times, tf_freqs, ersp_2, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('normalized - condition 2')

subplot(2, 3, 6)
contourf(tf_times, tf_freqs, ersp_1 - ersp_2, 40, 'linecolor','none')
clim([-6, 6])
colorbar()
title('normalized - condition difference')

colormap('parula')