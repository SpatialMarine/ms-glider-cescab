%% Conversion of spectra in hybrid millidecades
%
% Description:
% This function converts narrowband spectra in hybrid millidecades spectra. Those spectra have a
%   a frequency axis with a 1 Hz step up to 434 Hz with an increasing size afterwards. It is
%   compatible with 1D or 2D matrices of spectral values
%
% References:
% Martin, S. B. et al. (2021). Hybrid millidecade spectra: A practical format for exchange of
%   long-term ambient sound data, JASA Express Letters 1(1), 011203.
%   https://asa.scitation.org/doi/10.1121/10.0003324
% Martin, S. B. et al. (2021). Erratum: Hybrid millidecade spectra: A practical format for exchange
%   of long-term ambient sound data [JASA Express Letters 1(1), 011203 (2021)], JASA Express Letters
%   1(8), 081201. https://doi.org/10.1121/10.0005818
%
% Input arguments:
% - freq_nb          = Narrowband frequency axis, must have a frequency step of 1 Hz [Hz]
% - spectra_nb       = Narrowband spectral values [dB re 1μPa²/Hz]
%
% Output arguments:
% - freq_mdh         = Hybrid millidecades frequency axis [Hz]
% - spectra_mdh      = Hybrid millidecades spectral values [dB re 1μPa²/Hz]
% - limits_bands_hmd = Low and high limits of the HMD frequency bands [Hz]
%
% History:
% - 2024-05-24: Creation of the function (Pierre Mercure-Boissonnault)
% - 2024-10-17: Using the frequency axis of the erratum. (PMB)
% - 2026-06-04: Translation. (PMB)

function [freq_hmd, spectra_hmd, limits_bands_hmd] = ...
  convert_to_hybrid_millidecades(freq_nb, spectra_nb)
freq_threshold = 435.5; % [Hz]

% Check if the input axis have a step of 1 Hz
if ~all(diff(freq_nb)==1)
  error('The input frequency axis doesn''t have a constant frequency step of 1 Hz.');
end

% Selection of the narrowbands under the threshold
sel_freq_nb = freq_nb < freq_threshold;
freq_1Hz = freq_nb(sel_freq_nb);
spectra_1Hz = spectra_nb(sel_freq_nb, :);

limits_bandes_1Hz = [max(freq_1Hz-0.5,0) freq_1Hz+0.5];

% Computation of the millidecade part
% This uses the average_on_third_octaves function but with custom bands
[freq_md, spectra_md, limits_bands_md] = average_on_third_octaves(freq_nb, spectra_nb, ...
  false, 10^(1/1000), 1000/10^(1/2000), freq_threshold);

% Combining the two
freq_hmd = [freq_1Hz; freq_md];
spectra_hmd = [spectra_1Hz; spectra_md];
limits_bands_hmd = [limits_bandes_1Hz; limits_bands_md];
end
