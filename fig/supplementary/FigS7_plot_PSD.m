
%% Plot the power spectral densities of the CESCAB timeseries
%
% Description:
% This prepares the figures of the distribution of the noise recorded by the CESCAB mission.
%
% History:
% - 2025-10-23: Creation of the function (Greta Jankauskaite)
% - 2025-10-29: Optimization and argument for the histogram bins. (Pierre Mercure-Boissonnault)
% - 2026-06-04: Clean-up. (PMB)


% Give visibility to the subfunctions
addpath(genpath('Toolboxes'));

output_dir = fullfile('output', 'noise', 'fig');

% Create the output folder if it doesn't exist
if ~exist(output_dir, 'dir')
  mkdir(output_dir);
end

% Option to show the PSD plot in third-octave or hybrid millidecade bands
plot_third_octaves = true;

if plot_third_octaves % Third-octave values
  timeseries_data = load(fullfile('output', 'noise', 'timeseries', 'CESCAB_timeseries_to.mat'));
  freq_axis = timeseries_data.freq_to;
else % Hybrid millidecade values
  timeseries_data = load(fullfile('output', 'noise', 'timeseries', 'CESCAB_timeseries_hmd.mat'));
  freq_axis = timeseries_data.freq_hmd;
end

% Create a vector indicating what data to keep from the timeseries
sel_time = true(size(timeseries_data.time_axis));

% Recordings containing surfacing or glider self-noise are flagged in their filenames.
% A recording is retained only if its complete filename follows:
% CESCAB_yyyymmdd_HHMMSS.mat
%
% Files containing an additional suffix are excluded, for example:
%   CESCAB_20240801_121022_surface.mat
%   CESCAB_20240801_124022_glider_noise.mat
%   CESCAB_20240801_145522_surface_glider_noise.mat

file_list = dir(fullfile('output', 'MAT_30s', '*.mat'));
file_names = string({file_list.name});

is_clean = ~cellfun( ...
  'isempty', ...
  regexp(cellstr(file_names), '^CESCAB_\d{8}_\d{6}\.mat$', 'once') ...
);

noisy_file_names = file_names(~is_clean);

% Exclude the 30 s time interval associated with each noisy recording
for i_file = 1:length(noisy_file_names)
  file_start = datenum( ...
    extractBetween(noisy_file_names(i_file), 8, 22), ...
    'yyyymmdd_HHMMSS' ...
  );

  file_end = file_start + 30/(24*60*60);

  sel_time_file = timeseries_data.time_axis >= file_start & ...
                  timeseries_data.time_axis < file_end;

  sel_time(sel_time_file) = false;
end

%% Prepare the figure
fig = figure('WindowState', 'Maximized', 'color', 'w');
if plot_third_octaves
  psd_plot_function(freq_axis, timeseries_data.spectro_to(:, sel_time), fig, [60 130]);
else
  psd_plot_function(freq_axis, timeseries_data.spectro_hmd(:, sel_time), fig, [30 100]);
end


% Save the figure
if plot_third_octaves
  export_fig(fullfile(output_dir, 'CESCAB_PSD_third_octaves.png'), ...
    '-png', '-p0.01', '-m2', fig);
else
  export_fig(fullfile(output_dir, 'CESCAB_PSD_hybrid_millidecades.png'), ...
    '-png', '-p0.01', '-m2', fig);
end


% Local function to plot the power spectral densities and quantiles of a timeseries.
% Input arguments:
% - fHz       = Frequency axis [Hz]
% - dB        = Spectral values of the timeseries [dB]
% - fig       = Figure handle
% - dB_limits = Vector of the inferior and superior limits of the histogram bins [dB].
%
function psd_plot_function(fHz, dB, fig, dB_limits)
% Prepare the power spectral density values
dBstep = .1;
dBscale = dB_limits(1):dBstep:dB_limits(2); % Prepare the bins
for ifr = 1:length(fHz)
  psd_values(ifr, :) = histcounts(dB(ifr,:), [dBscale dBscale(end)+dBstep], 'Normalization', 'pdf');
end
psd_values(psd_values==0) = NaN;

figure(fig);

% Plot the psd
pcolor(fHz, dBscale, psd_values');
shading flat;

% Prepare the colorbar
cb = colorbar;
cb.Label.String = 'Empirical Probability Density';
caxis([0 0.1]);

% Format the frequency axis
xlim([10 70000]);
set(gca, 'xscale', 'log', 'TickDir', 'out');
exponents_ticks = (ceil(log10(min(xlim))*3):floor(log10(max(xlim))*3))/3;
xticks(round(10.^exponents_ticks, 1, 'significant'));

% Compute the quantiles
quantile_axis = [.05 .25 .5 .75 .95];
quantiles_curves = quantile(dB, quantile_axis, 2);

% Plot the quantile curves
hold on;
for i_q = 1:length(quantile_axis)
  plot(fHz,quantiles_curves(:,i_q),"linewidth",2);
end
hold off;

% Additional formatting
xlabel('Frequency (Hz)')
ylabel('PSD (dB re 1 $\mu$Pa$^2$ Hz$^{-1}$)',"Interpreter","latex")
grid on;
box on;
end