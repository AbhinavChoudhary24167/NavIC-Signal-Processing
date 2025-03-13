clear;
close all;
clc;

%% System Parameters
c_len = 1023; % CA code length
satellite_no = 7; % Number of NavIC satellites
fs_tx = 1.023 * 1e6; % Sampling rate at transmitter, 1.023 MHz
fs_rx = 2 * 1e6; % Sampling rate at receiver, 2 MHz 
T = 1e-3; % Time duration
c = 3e8; % Speed of light
Nsamp_nu = round(fs_rx * T); % Total Number of Samples
fd = -6e3 : 500 : 6e3; % Doppler Frequency axis
delt = 1/fs_rx; % Sampling time
time_axis = (0:Nsamp_nu-1) * delt; % Time axis (seconds)
range_axis = c * time_axis; % Range axis (meters)
nfft = 2048; % FFT size
detection_threshold_nu = 3; % Peak detection threshold

%% Load the NavIC PRN codes of all satellites
navic_prn = load('navic_prn.mat').navic_prn;

% Resample PRN sequence at receiver
fracRepFactor = fs_rx / fs_tx;
[upfac, downfac] = rat(fracRepFactor, 1e-8 * norm(fracRepFactor,1));
idx = 1:downfac:upfac*size(navic_prn,1);
ca_code = navic_prn(ceil(idx/upfac), :); % Resampled PRN codes

%% Load GNSS Simulator datasets
N = 5; % Number of test cases
for j = N
    disp(['\n=== Test Case:', num2str(j), ' ===']);
    file = strcat('in_data', string(j), '.mat');
    in = load(file).in;

    % Construct Doppler phase shifts
    phase = exp((-1 * 2 * pi * 1j * time_axis).' .* fd);

    % FFT of transmit PRN sequence
    prn_zp = zeros(nfft, satellite_no);
    prn_zp(1:Nsamp_nu, :) = ca_code; % Zero padding to nfft
    fft_prn = conj(fft(prn_zp, nfft)); % NFFT-point FFT for each satellite  

    % Signal Acquisition for each satellite
    sel_satind = 0;
    max_peak_mag = -Inf;  % Track the strongest detected satellite
    
    for prnid = 1:satellite_no

        % RSP function call for range-Doppler ambiguity
        [mag_out_rd(:,:,prnid), peak_mag, peak_idx] = RSP(fd, phase, nfft, Nsamp_nu, fft_prn(:,prnid), in);
        
        % Debugging: Print peak magnitude from RSP
        disp(['Satellite ', num2str(prnid), ' Peak Mag from RSP: ', num2str(peak_mag)]);
        
        % Threshold comparison
        y = mag_out_rd(:,peak_idx,prnid); % Extract the slice
        [peaks, locs] = findpeaks(y, (1:nfft)); % Finding all local maxima
        
        % Find the index of the maximum peak within the peaks array
        idx_peak = find(peaks == peak_mag, 1);
        
        % Print additional peak details (if available)
        if isempty(idx_peak)
            disp(['Satellite ', num2str(prnid), ': idx_peak is empty.']);
        elseif idx_peak <= 1 || idx_peak >= length(peaks)
            disp(['Satellite ', num2str(prnid), ': idx_peak at boundary, cannot print adjacent peaks.']);
        else
            disp(['Satellite ', num2str(prnid), ' Details:']);
            disp(['    peak_mag: ', num2str(peak_mag)]);
            disp(['    Previous Peak (peaks(idx_peak-1)): ', num2str(peaks(idx_peak-1))]);
            disp(['    Next Peak (peaks(idx_peak+1)): ', num2str(peaks(idx_peak+1))]);
        end
        
        if isempty(peaks)
            disp(['Satellite ', num2str(prnid), ': No peaks found! Skipping...']);
            continue; % Move to the next satellite
        end
        
        % Ensure idx_peak is valid for detection threshold check
        if isempty(idx_peak) || idx_peak <= 1 || idx_peak >= length(peaks)
            detected = false;
        else
            detected = ((peak_mag / peaks(idx_peak+1) > detection_threshold_nu) || ...
                        (peak_mag / peaks(idx_peak-1) > detection_threshold_nu));
        end

        if detected
            disp(['Satellite ', num2str(prnid), ': Detection Threshold Met']);
            if peak_mag > max_peak_mag
                max_peak_mag = peak_mag;
                sel_satind = prnid;
            end
        else
            disp(['Satellite ', num2str(prnid), ': Detection Not Met']);
        end
    end
writematrix(mag_out_rd,"mag_is_OUTtoo_2048.csv");

    % Finding range and Doppler of peak for acquired satellite
    if sel_satind > 0
        out = mag_out_rd(:,:,sel_satind);
        out_peak = max(max(out));
        [range_idx, doppler_idx] = find(out == out_peak);

        % Print detected range and Doppler frequency
        detected_range = range_axis(range_idx);
        detected_doppler = fd(doppler_idx);
        disp(['✅ Test Case ', num2str(j), ':']);
        disp(['  - Detected Satellite: ', num2str(sel_satind)]);
        disp(['  - Detected Range: ', num2str(detected_range), ' meters']);
        disp(['  - Detected Doppler Frequency: ', num2str(detected_doppler), ' Hz']);
    else
        disp(['❌ Test Case ', num2str(j), ': No satellite detected.']);
    end

    % Plotting range-Doppler image of detected satellite
    if sel_satind > 0
        figure;
        imagesc(fd, range_axis, out(1:Nsamp_nu,:), [max(max(out))-20, max(max(out))]);
        colorbar;
        xlabel('Doppler Frequency (Hz)');
        ylabel('Range (m)');
        title(['Satellite ', num2str(sel_satind), ' Range-Doppler Plot']);
        colormap(jet(256));
        set(gca, 'YDir', 'normal');
    end
end