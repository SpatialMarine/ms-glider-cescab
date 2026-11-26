
%% Preparation of timeseries for the CESCAB recordings
% Description:
% This script combines the spectral values previously calculated for each
% 30 s recording of the CESCAB mission into a third-octave-band timeseries.
%
% History:
% - 2025-08-28: Creation of the script (Pierre Mercure-Boissonnault)
% - 2026-06-04: Clean-up and translation (PMB)

input_dir = fullfile('output', 'MAT_30s');
output_dir = fullfile('output', 'noise', 'timeseries');

list_mat_files = dir(fullfile(input_dir, '*.mat'));
nb_mat_files = length(list_mat_files);

% Create the output folder if it doesn't exist
if ~exist(output_dir, 'dir')
  mkdir(output_dir);
end

%% Load and combine the third-octave data
spectro_to = [];
time_axis = [];

for i_file = 1:nb_mat_files
  % Load the file
  spectro_data = load(fullfile( ...
    list_mat_files(i_file).folder, ...
    list_mat_files(i_file).name));

  % Append the data at the end of the matrices
  new_indexes_time = (1:length(spectro_data.time_axis)) + length(time_axis);
  time_axis(new_indexes_time) = spectro_data.time_axis;
  spectro_to(:, new_indexes_time) = spectro_data.spectro_to;
end

%% Save the third-octave timeseries
freq_to = spectro_data.freq_to;

save(fullfile(output_dir, 'CESCAB_timeseries_to.mat'), ...
  'spectro_to', 'freq_to', 'time_axis');