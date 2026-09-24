%% Average of narrowband spectral values on third octave bands
% Description:
% This function averages the narrowbands values on third octaves by integrating over each band. The
%   function also allows to do the operation on decidecade bands or on custom (using base_type and
%   freq_ref).
%
% Input arguments:
% - freq_nb        = Narrowband frequency axis [Hz]
% - spectra_nb     = Narrowband spectral values [dB re 1μPa²/Hz]
% - value_is_in_dB = Boolean indicating if the spectra values are in dB (requiring a conversion in
%                      power units, optional, true by default)
% - base_type      = Choice of the factor between 2 bands (optional)
%                    - false           : base 2 (by default)
%                    - true            : base 10 (decidecades)
%                    - numerical value : custom factor between the bands
% - freq_ref       = Reference frequency (optional, 1000Hz by default) [Hz]
% - freq_min       = Minimum frequency to keep (optional, starts with the first non-zero
%                      frequency) [Hz]
%
% Output arguments:
% - freq_to         = Third octave frequency axis [Hz]
% - spectra_to      = Third octave spectral values [dB re 1μPa]
% - limits_bands_to = Low and high limits of the third octave bands [Hz]
%
% History:
% - 2022-04-07: Creation of the function (Pierre Mercure-Boissonnault)
% - 2022-08-30: Possibility to choose custom bands. (PMB)
% - 2024-05-24: Added the argument freq_min. (PMB)
% - 2025-04-15: Optimisation for large matrices. (PMB)
% - 2026-06-04: Clean-up and translation. (PMB)

function [freq_to, spectra_to, limits_bands_to] = ...
  average_on_third_octaves(freq_nb, spectra_nb, value_is_in_dB, base_type, freq_ref, freq_min)
% Preparing default values for optionnal arguments
if ~exist('value_is_in_dB', 'var')
  value_is_in_dB = true;
end
if ~exist('base_type', 'var')
  base_type = [];
end
if ~exist('freq_ref', 'var')
  freq_ref = [];
end
if ~exist('freq_min', 'var')
  freq_min = min(freq_nb(freq_nb>0));
end
  
% Preparing the third octave frequency axis
freq_range = [freq_min max(freq_nb)];
[freq_to, width_bands, limits_bands_to] = ...
  prepare_third_octave_axis(freq_range, base_type, freq_ref);

% Rearranging the dimensions to make sure that it starts with the frequencies
sz = size(spectra_nb);
if length(sz) == 2 && sz(1) == 1 && sz(2) > 1
  spectra_nb = spectra_nb';
  sz = flip(sz);
end

% Conversion of the spectral values in power if needed
if value_is_in_dB
  spectra_nb = 10.^(spectra_nb/10);
end

% Interpolation of the narrowband spectra over the third octave band limits
% Note : It is looking for the frequencies before each limit to do the interpolation only on the
%   value before and after each limit (instead of doing the interpolation using the whole set of
%   values)
limits_bands_vector = [limits_bands_to(:,1); limits_bands_to(end,2)];
i_freq_prec_limits = interp1(freq_nb, 1:length(freq_nb), limits_bands_vector, 'previous');
i_freq_list_interp = reshape(i_freq_prec_limits' + [0; 1], 2*length(i_freq_prec_limits), 1);
i_freq_list_interp = unique(i_freq_list_interp); % Remove duplicates
values_limits_bands = interp1(freq_nb(i_freq_list_interp), spectra_nb(i_freq_list_interp, :), ...
  limits_bands_vector);

spectra_to = NaN([length(freq_to) sz(2:end)]);

% Integration over each band
for i_band = 1:length(freq_to)
  % Prepare the frequency axis for the integration
  sel = (i_freq_prec_limits(i_band)+1):i_freq_prec_limits(i_band+1);
  freq_integration = [limits_bands_to(i_band, 1) freq_nb(sel)' limits_bands_to(i_band, 2)];

  % Prepare the list of values to integrate
  valeurs_integration = [values_limits_bands(i_band, :); spectra_nb(sel, :); ...
    values_limits_bands(i_band+1, :)];

  % Integration
  spectra_to(i_band,:) = trapz(freq_integration, valeurs_integration(:,:))/width_bands(i_band);
end

% Reconversion in dB if needed
if value_is_in_dB
  spectra_to = 10*log10(spectra_to);
end
end
