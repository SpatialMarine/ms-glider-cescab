
%% Preparation of the LTSA spectrograms
% Description:
% This script prepares figures of the LTSA spectrograms of the 30s recordings of the CESCAB mission.
%
% History:
% - 2026-01-14: Creation of the script (Pierre Mercure-Boissonnault)
% - 2026-06-04: Clean-up. (PMB)

% Add the paths the subfunctions and toolboxes
addpath(genpath('Subfunctions'));
addpath(genpath('Toolboxes'));


input_dir = fullfile('output', 'MAT_30s');
output_dir = fullfile('output', 'noise', 'fig');

% Create the output folder if it doesn't exist
if ~exist(output_dir, 'dir')
  mkdir(output_dir);
end

% Options
% Target number of points/pixels to plot on the time axis
target_time_points = 2000;

% Format the frequency axis based on the progression of the hybrid millidecade scale (linear, then
%   log). If false, use a log scale.
freq_scale_hmd = true;

% Start and end dates
% date_range = [-Inf Inf]; % Use to whole dataset
date_range = [datenum([2024 08 08 18 0 0]) datenum([2024 08 09 6 0 0])]; % 12h period

duration = diff(date_range); % [days]
est_number_of_recordings = duration*24*60/7.5;
block_size_ini = max(round(59*est_number_of_recordings/target_time_points), 1);

%% Loading the data
file_list = dir(fullfile(input_dir, '*.mat'));

% Select only the recordings within the data range
start_dates_30s = datenum(extractBetween({file_list.name}, 8, 22), 'yyyymmdd_HHMMSS');
file_list = file_list(start_dates_30s >= date_range(1) & start_dates_30s < date_range(2));


% Recordings containing surfacing or glider self-noise are flagged in their filenames.
% A recording is retained only if its complete filename follows:
% CESCAB_yyyymmdd_HHMMSS.mat
%
% Files containing an additional suffix are excluded, for example:
%   CESCAB_20240801_121022_surface.mat
%   CESCAB_20240801_124022_glider_noise.mat
%   CESCAB_20240801_145522_surface_glider_noise.mat

file_names = string({file_list.name});

is_clean = ~cellfun( ...
  'isempty', ...
  regexp(cellstr(file_names), '^CESCAB_\d{8}_\d{6}\.mat$', 'once') ...
);

is_noisy = ~is_clean;

time_axis = [];
time_axis_plot = [];
ltsa_spectrogram = [];
i_pt = 0;
i_pt_plot = 0;

tic;
for i_file = 1:length(file_list)
  % Reading the .mat spectrogram of each file
  spectrogram_data = load([file_list(i_file).folder filesep file_list(i_file).name]);
  if i_file == 1
    freq_axis = spectrogram_data.freq_hmd;
  end

  % Determining the number of blocks to split the recording
  nb_windows = length(spectrogram_data.time_axis);
  nb_blocks = max(ceil(nb_windows/block_size_ini), 1);

  if ~is_noisy(i_file) % If the file is clean 
    % Readjustment of the block size to reduce the possibility of having smaller blocks at the end
    block_size = ceil(nb_windows/nb_blocks);

    for i_block = 1:nb_blocks
      % Preparation of the indexes to select the data of the block
      i_sel = (1+(i_block-1)*block_size):min((i_block*block_size), nb_windows);

      % Averaging the time and spectral values on the block
      i_pt = i_pt + 1;
      i_pt_plot = i_pt_plot + 1;
      time_axis(i_pt) = mean(spectrogram_data.time_axis(i_sel));
      time_axis_plot(i_pt) = i_pt_plot;
      ltsa_spectrogram(:, i_pt) = ...
        10*log10(mean(10.^(spectrogram_data.spectro_hmd(:, i_sel)/10), 2));
    end
  else
    % Skip the noisy file while preserving its space on the time axis
    i_pt_plot = i_pt_plot + nb_blocks;
  end
end


toc; % It's taking me ~5min on the whole dataset, it should be faster for you as I'm slowed down by
%   reading the files on the network

tic;

fig = figure('Color', 'w', 'WindowState', 'Maximized');% Create a new figure in full screen

% Set the aspect ratio of the axes if needed
% axis('square');
% set(gca, 'PlotBoxAspectRatio', [1 3/4 1], 'DataAspectRatioMode', 'auto');

% Adding an extra point on the time axis so that the last spectrogram column will show in the pcolor
time_axis(end+1) = time_axis(end) + mean(diff(time_axis));
time_axis_plot(end+1) = time_axis_plot(end) + 1;
ltsa_spectrogram(:, end+1) = NaN(size(ltsa_spectrogram, 1), 1);

% Adding an extra frequency so that the last spectrogram row will show in the pcolor
freq_axis(end+1) = freq_axis(end)*10^(-1/1000);
ltsa_spectrogram(end+1,:) = NaN(1, size(ltsa_spectrogram, 2));

% Setting the frequency axis used in the pcolor
if freq_scale_hmd
  freq_axis_plot = 1:length(freq_axis);
else
  % Preparing the lower limits of the frequency axis as it is was pcolor is using in 'flat' mode
  lower_freq_limit = freq_axis - 0.5;
  lower_freq_limit(1) = 0;
  sel_log_part = [false; diff(freq_axis)>1];
  lower_freq_limit(sel_log_part) = freq_axis(sel_log_part)*10^(-1/2000);
  freq_axis_plot = lower_freq_limit/1000;
end

% Plotting the spectrogram
pcolor(time_axis_plot, freq_axis_plot, ltsa_spectrogram);
shading flat;

% Setting the x and y labels
duration = diff(time_axis([1 end])); % Updating the duration
if duration < 1
  xlabel('Time');
else
  xlabel('Date');
end
ylabel('Frequency [kHz]');

% Preparing the colorbar
cmap = imread('colormap plasma raven.png');
cmap = double(flip(permute(cmap, [1 3 2])))/255; % Permutation, inversion and normalization
colormap(cmap);
caxis([50 80]);
cb = colorbar();
cb.Label.String = 'Received Levels [dB re 1μPa²/Hz]';

grid on;
box on;
set(gca, 'Layer', 'top');
set(gca, 'FontSize', 13);

xlim(time_axis_plot([1 end]));

% Preparation of the frequency axis
if freq_scale_hmd
  ylim(freq_axis_plot([1 end]));

  % Setting the ticks
  freq_ticks = [0 0.2 0.5 1 2 5 10 20 50 100];
  sel_log_part = freq_ticks > 0.435; % or 0.455 depending on how the hmd are defined
  i_freq_ticks = [interp1(freq_axis/1000, freq_axis_plot+0.5, freq_ticks(~sel_log_part)) ...
    interp1(log(freq_axis(2:end)/1000), freq_axis_plot(2:end)+0.5, log(freq_ticks(sel_log_part)))];
  yticks(i_freq_ticks);
  yticklabels(cellstr(string(freq_ticks)));
else
  ylim([0.01 freq_axis(end)/1000]);
  set('Yscale', 'log');
  marques_frequence('y');
end

% Preparation of the time axis ticks
if duration < 0.5
  ticks_by_day = 24;
  tick_format = 'HH:MM';
elseif duration < 1
  ticks_by_day = 12;
  tick_format = 'HH:MM';
elseif duration < 2
  ticks_by_day = 4;
  tick_format = 'dd mmm HH:SS';
elseif duration < 2.5
  ticks_by_day = 3;
  tick_format = 'dd mmm HH:SS';
elseif duration < 3
  ticks_by_day = 2;
  tick_format = 'dd mmm HH:SS';
else
  ticks_by_day = 1/max(round(duration/10));
  tick_format = 'dd mmm';
end

time_ticks = (ceil(time_axis(1)*ticks_by_day):floor(time_axis(end)*ticks_by_day))/ticks_by_day;
i_time_ticks = interp1(time_axis, time_axis_plot+0.5, time_ticks);
xticks(i_time_ticks);
xticklabels(datestr(time_ticks, tick_format));
toc; % <10s

% Save the figure with magnification
export_fig(fullfile(output_dir, 'LTSA_spectrogram.png'), ...
  '-png', '-p0.01', '-m2', fig);