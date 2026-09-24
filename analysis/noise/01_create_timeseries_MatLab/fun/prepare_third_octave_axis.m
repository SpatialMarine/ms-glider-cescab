%% Preparation of the third octave axis
% Description:
% This function prepares a third octave frequency axis. It can also prepare a decidecade axis or a
%   custom one.
%
% References:
% ASA/ANSI S1.11-2004 (R2009), American National Standard Specification for Octave-Band and
%   Fractional-Octave-Band Analog and Digital Filters
%
% Input arguments:
% - freq_range = Limits of the frequency axis wanted ([fmin fmax]) [Hz]. The frequencies must be
%                  positive and non-zero. The third octave bands must be fully contained in the
%                  frequency range to be returned.
% - base_type  = Choice of the factor between 2 bands (optional)
%                - false           : base 2 (by default)
%                - true            : base 10 (decidecades)
%                - numerical value : custom factor between the bands
% - freq_ref   = Reference frequency (optional, 1000Hz by default) [Hz]
%
% Output arguments:
% - freq_to      = Third octave frequency axis [Hz]
% - width_bands  = Width of the third octave bands [Hz]
% - limits_bands = Low and high limits of the third octave bands [Hz]
%
% History:
% - 2022-02-16: Creation of the function (Pierre Mercure-Boissonnault)
% - 2022-04-13: Argument to choose the limits of the frequency axis. (PMB)
% - 2022-08-30: Added the possibility to prepare custom bands. (PMB)
% - 2026-06-04: Clean-up and translation. (PMB)

function [freq_to, width_bands, limits_bands, factor_bands] = ...
  prepare_third_octave_axis(freq_range, base_type, freq_ref)
% Checking the frequency axis
if any(freq_range <= 0)
  error('Invalid frequencies range values')
end

% Prepare default values for optionnal arguments
if ~exist('base_type', 'var') || isempty(base_type)
  base_type = false; % Base 2 by default
end
if ~exist('freq_ref', 'var') || isempty(freq_ref)
  freq_ref = 1000; % 1000Hz by default (used for both base 2 and 10)
end

% Set the factor between each band
if islogical(base_type)
  if ~base_type
    factor_bands = 2^(1/3); % Base 2
  else
    factor_bands = 10^(1/10); % Base 10
  end
else
  factor_bands = base_type; % Custom factor between bands
end

% Prepare the list of indexes for the bands (based on the ISO standard)
% Note: Only takes the bands that are fully covered by the frequency range
index_ini = log(freq_range(1)/freq_ref)/log(factor_bands)+30;
index_ini = ceil(index_ini + 0.5);
index_end = log(freq_range(2)/freq_ref)/log(factor_bands)+30;
index_end = floor(index_end - 0.5);
indexes_bands = (index_ini:index_end)';

% Prepare the frequency axis, limits and band widths
freq_to = freq_ref*factor_bands.^(indexes_bands-30);
limits_bands = freq_to.*[1/sqrt(factor_bands) sqrt(factor_bands)];
width_bands = diff(limits_bands, 1, 2);
end

