
# Liao. et al. 2022 Geoscientific Model Development

**Representing lateral groundwater flow in Earth system models**

Chang Liao<sup>1\*</sup>,
L. Ruby Leung<sup>1</sup>,
Yilin Fang<sup>2</sup>,
Teklu Tesfa<sup>2</sup>,
Robinson Negron-Juarez<sup>3</sup>,


<sup>1 </sup> Atmospheric Sciences and Global Change, Pacific Northwest National Laboratory, Richland, WA, USA

<sup>2 </sup> Hydrology Group, Pacific Northwest National Laboratory, Richland, WA, USA

<sup>3 </sup> Climate and Ecosystem Sciences Division, Lawrence Berkeley National Laboratory, Berkeley, CA, USA



\* corresponding author:  chang.liao@pnnl.gov

## Abstract

Lateral groundwater flow plays an important role in controlling water table dynamics. Due to the relatively coarse spatial resolutions of land surface models, this process is often omitted even though it can be significant due to subgrid heterogeneity. In this study, we developed a physical based model to simulate lateral groundwater flow using hillslopes to represent subgrid spatial variability in topography. This model explicitly considers the smooth transition between different water table scenarios (e.g., with or without a seepage phase in the lower elevations). We coupled this model to the land component (ELM) and river component (MOSART) of the Energy Exascale Earth System Model (E3SM) and applied it at global scale. Simulations show that lateral groundwater flow is affected by both topography and river stage through their impacts on water table gradient. Future improvements are needed because representation of hillslope is the key to model lateral flow processes.

## Journal reference
Liao. et al. (2022). Representing lateral groundwater flow in Earth system models.

## Contributing modeling software

| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| HexWatershed | 3.0 | https://github.com/changliao1025/pyhexwatershed | doi.org/10.5281/zenodo.6425881 |
| PyEarth | 0.1.25 | https://github.com/changliao1025/pyearth | doi.org/10.5281/zenodo.6368652 |
| PyE3SM | 0.1.0 | https://github.com/changliao1025/pye3sm | doi.org/10.5281/zenodo.7591982 |
| E3SM | 0.1.0 | https://github.com/changliao1025/E3SM/tree/changliao/elm/hillslope |  |

## Reproduce my experiment

You need to follow three major steps to reproduce this study:

1. Install the pyhexwatershed software
2. Download the DEM for the study area
3. Run the HexWatershed model to generate the hillslope data
4. Install the E3SM model with the hillslope capability
5. Install the PyE3SM software
6. Run the PyE3SM to setup the E3SM model with the hillslope capability
7. Run the E3SM model

## Reproduce my figures

Use the scripts found in the `codes/k34/analysis/plot/` directory to reproduce the figures used in this publication.


