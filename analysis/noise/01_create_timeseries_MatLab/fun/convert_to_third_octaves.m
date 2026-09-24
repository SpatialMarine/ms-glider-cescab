%% Conversion of narrowband spectra in third octave
% Description:
% This function converts narrowband spectra in third octave or decidecade spectra. It is compatible
%   with 1D or 2D matrices of spectral values.
%
% Input arguments:
% - freq_nb    = Narrowband frequency axis [Hz]
% - spectra_nb = Narrowband spectral values [dB re 1μPa²/Hz]
% - base_type  = Choice of the factor between 2 bands (optional)
%                - false           : base 2 (by default)
%                - true            : base 10 (decidecades)
%                - numerical value : custom factor between the bands
% - freq_ref   = Reference frequency (optional, 1000Hz by default) [Hz]
% - freq_min   = Minimum frequency to keep [Hz]
%
% Output arguments:
% - freq_to         = Third octave frequency axis [Hz]
% - spectra_to      = Third octave spectral values [dB re 1μPa]
% - limits_bands_to = Low and high limits of the third octave bands [Hz]
%
% History:
% - 2022-04-07: Creation of the function (Pierre Mercure-Boissonnault)
% - 2023-01-27: Simplifying the function using average_on_third_octaves. (PMB)
% - 2026-06-04: Translation. (PMB)

function [freq_to, spectra_to, limits_bands_to] = ...
  convert_to_third_octaves(freq_nb, spectra_nb, varargin)
  
% Compute the average of the third octave bands
[freq_to, spectra_to, limits_bands_to] = ...
  average_on_third_octaves(freq_nb, spectra_nb, false, varargin{:});

% Conversion into a sum to get power values
width_bands = diff(limits_bands_to, 1, 2);
for i_band = 1:length(freq_to)
  spectra_to(i_band,:) = spectra_to(i_band,:) + 10*log10(width_bands(i_band));
end
end
