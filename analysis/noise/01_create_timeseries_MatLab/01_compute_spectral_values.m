
%% Computation of the spectral values of the CESCAB recordings
% Description:
% This script prepares the spectral values for each 30s recording of the CESCAB mission in third
%   octave bands and in hybrid millidecade bands.
%
% History:
% - 2025-08-28: Creation of the script (Pierre Mercure-Boissonnault)
% - 2026-06-04: Clean-up and translation. (PMB)

% path to custom functions
addpath(genpath(fullfile('analysis', 'noise', '01_matlab', 'fun')));

tic_tot = tic; % Start a timer for the overall time taken by the script

input_dir = fullfile('input', '1_CESCAB-raw-flac');
output_dir = fullfile('output', 'MAT_30s');

% List the FLAC files
list_flac = dir(fullfile(input_dir, '**', '*.flac'));
list_flac = list_flac(~startsWith({list_flac.name}, '.'));

% Create the output folder if it doesn't exist
if ~exist(output_dir, 'dir')
  mkdir(output_dir);
end

% Loop over each recording
for i_file = 1:length(list_flac)
  % Prepare the path for the output file
  filename = list_flac(i_file).name;
  out_path = fullfile(output_dir, [filename(1:end-5) '.mat']);

  % Read the recording
  [signal, sampling_rate] = audioread([list_flac(i_file).folder filesep filename]);

  % Prepare the narrowband psd values of the spectrogram 
  nfft = sampling_rate;
  fft_window = hann(nfft);
  overlap = nfft/2;
  [~, freq_axis_nb, time_axis, fft_amplitudes] = ...
    spectrogram(signal, fft_window, overlap, nfft, sampling_rate, 'psd');

  % Convert the result of the FFT in dB
  hydrophone_sensitivity = -206;
  V_peak = 1.5;
  gain = 13.2 + 20; % Including preamp gain
  sensitivity_end_to_end = -20*log10(V_peak) + hydrophone_sensitivity + gain;
  spectro_nb = 10*log10(fft_amplitudes) - sensitivity_end_to_end;

  % Convert the narrowbands in hybrid millidecades
  [freq_hmd, spectro_hmd] = convert_to_hybrid_millidecades(freq_axis_nb, spectro_nb);

  % Convert the narrowbands in third octave (starting from 10 Hz)
  [freq_to, spectro_to] = convert_to_third_octaves(freq_axis_nb, spectro_nb, false, 1000, 8.5);

  % Prepare additionnal time information
  start_date = datenum(filename(8:22), 'yyyymmdd_HHMMSS');
  duration = length(signal)/sampling_rate;
  time_axis = time_axis/24/3600 + start_date;

  % Save the spectal values in a .mat file
  save(out_path, 'spectro_hmd', 'freq_hmd', 'spectro_to', 'freq_to', 'time_axis', 'duration', ...
    'start_date');
end
disp(['Total time used : ' num2str(toc(tic_tot)) 's']);