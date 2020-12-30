function [omegas, M_psd] = power_spectrum(Mitral_spike_history, dt, num_averages) 

% SS1 = squeeze(Mitral_spike_history(:, 1, ceil(end/2):end));
SS1 = Mitral_spike_history(:, ceil(end/2):end);

% number of chunks the spike data is broken into
num_steps = size(SS1, 2);

Nm = size(SS1, 1);

% size of those chunks
% the larger this is, the higher resolution you get at low freqs
num_fft_steps = floor(num_steps/num_averages);

% occasionally this bullshit is necessary
if mod(num_fft_steps, 2)==1
    num_fft_steps = num_fft_steps-1;
end

% unsure of significance of sorting but afraid to take it out 
[~, sort_indices] = sort(sum(SS1, 2), 'descend');
SS1_sorted = SS1(sort_indices, :);

% I renamed all these variables because their previous names
% didn't make sense to me.
M_psd = zeros(Nm, num_fft_steps);

% for each chunk
for fft_i = 1:num_averages
    % grab a time chunk of the data 
    nti_min = floor((fft_i-1)*num_steps/num_averages) + 1;
    nti_max = min(nti_min + num_fft_steps - 1, size(SS1_sorted, 2));
    fourier_chunk = SS1_sorted(:, nti_min:nti_max);

    % subtract the mean and fft it
    % Matlab's fft includes the 2pi term, so this is in radians
    M_fft = abs(fft(fourier_chunk-repmat(mean(fourier_chunk, 2), 1, size(fourier_chunk,2)), ...
        num_fft_steps, 2));

    % in calculating the PSD, normalize by 
    %   1) length of sample 
    %   2) sampling frequency (dt)
    % if using ms data, units are already in milliseconds
    M_psd = M_psd + dt*M_fft.^2/num_fft_steps;
end

% need to compute the power inside the loop and THEN take the avg
% otherwise noise terms get muted and power spectrum will be off
% by autocorrelation of the noise. see mrm_spectral.m
M_psd = M_psd/num_averages;

max_period = num_fft_steps*dt; % maximum period, in ms
min_freq = 1/max_period; % corresponding minimum frequency, in kHz
omegas = min_freq*[0:num_fft_steps/2, -num_fft_steps/2+1:-1];

end