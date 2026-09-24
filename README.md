## Ocean Gliders for Passive Acoustic Monitoring of Cetaceans and Underwater Noise in a large Marine Protected Area - R Code

This repository provides the R code used in the following research paper:

Greta Jankauskaite <a href="https://orcid.org/0000-0002-2177-6815"><img src="images/orcid.svg" alt="ORCID logo" width="21"/></a>, Simone Antichi, Pierre Cauchy, Margalida Cerdá, Marissa Garcia, Holger Klinck, Pierre Mercure-Boissonnault, Albert Miralles, Joaquín Tintoré, Thomas Webber, Nikolaos Zarokanellos, Juan Antonio Raga, David March - *Ocean gliders for passive acoustic monitoring of cetaceans and underwater noise in a large marine protected area*, doi:

------------------------------------------------------------------------

### Repository structure

| Folder | Description |
|----|----|
| `analysis/cetaceans/01_inter-annotator` | Inter-annotator comparison, creation of the final detection table. |
| `analysis/cetaceans/02_detection_statistics` | Calculate detection percentages and diel patterns analysis. |
| `analysis/noise` | Calculate sound pressure levels (SPL) in the 63 Hz and 125 Hz third-octave bands used to assess shipping noise, and the 8 kHz band used to assess wind-related noise. | |
| `setup.R` | Environment configuration, package dependencies, and global parameters. |

### Data availability

The acoustic data required to reproduce the analyses are publicly available on Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20702261.svg)](https://doi.org/10.5281/zenodo.20702261)

### License

Copyright (c) 2026. Greta Jankauskaite

@gretajan97 <a href="https://github.com/gretajan97"><img src="images/github-logo-black.svg" alt="GitHub logo" width="23"/></a> <a href="https://www.linkedin.com/in/greta-jankauskaite-510113119/"><img src="images/linkedin_short.svg" alt="Linkedin logo" width="23"/></a>

Licensed under the [MIT License](https://github.com/SpatialMarine/ms-glider-cescab/blob/main/LICENSE)

### Acknowledgements

This study was supported by the CESCAB project (SOCIB Glider Open Access, SOCIB 2024 N-1 AAC 0040) and TECMAR project (Organismo Autónomo Parques Nacionales, 3002/2023). GJ was supported by a predoctoral grant of the Conselleria de Innovación, Universidades, Ciencia y Sociedad Digital (Generalitat Valenciana) (CIACIF/2021/049), as well as mobility grants from Generalitat Valenciana (CIBEFP/2023/082,2024/153) and FRQNT PBEEE (371612) for international research and training visits. DM acknowledges support from the CIDEGENT program of the Generalitat Valenciana (CIDEGENT/2021/058). We thank crew, engineers, and technicians who participated in the ocean glider deployment and piloting.

