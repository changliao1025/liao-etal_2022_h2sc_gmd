[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14003483.svg)](https://doi.org/10.5281/zenodo.14003483)


# Liao. et al. 2024 Geoscientific Model Development

**Representing lateral groundwater flow between land and river in Earth system models**

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
Liao. et al. (2024). Representing lateral groundwater flow between land and river in Earth system models.

## Contributing modeling software

| Model | Version | Repository Link | DOI |
|-------|---------|-----------------|-----|
| HexWatershed | 3.0 | https://github.com/changliao1025/pyhexwatershed | doi.org/10.5281/zenodo.6425881 |
| PyEarth | 0.1.25 | https://github.com/changliao1025/pyearth | doi.org/10.5281/zenodo.6368652 |
| PyE3SM | 0.1.0 | https://github.com/changliao1025/pye3sm | doi.org/10.5281/zenodo.7591982 |
| E3SM | 0.1.0 | https://github.com/changliao1025/E3SM/tree/changliao/elm/hillslope |  |

## Reproduce my experiment

To reproduce the hillslope definition simulations, you need to follow the steps below:
1. Install the pyhexwatershed software. Please refer to HexWatershed documentation (https://hexwatershed.readthedocs.io/en/latest/) for more details.
2. Download the DEM for the study area. This is done using the Latitude = -2.6091 and Longitude = -60.2093 to define a 0.5 degree by 0.5 degree box. Then users can use the box to download the high resolution DEM from other resources such as the SRTM data. In this study, we use the Python package `elevation` (https://pypi.org/project/elevation/) to download the DEM data. An example box geojson file is provided in the `data/` directory.
3. Set up and run the HexWatershed model to generate the hillslope data. In thi case, the HexWatershed model is run in the elevation-only mode. See documentation for more details (https://hexwatershed.readthedocs.io/en/latest/). An example HexWatershed configuration file is provided in the `data/` directory.


To reproduce the E3SM simulations, you need to follow the steps below:
1. Install the E3SM model with the hillslope capability (https://github.com/changliao1025/E3SM/tree/changliao/elm/hillslope)
2. Install the PyE3SM package (https://github.com/changliao1025/pye3sm)
3. Install the PyEarth software (https://pypi.org/project/pyearth/)
4. Download the default E3SM ELM input data, including the forcing data such as precipation and temperature, and the surface/domain data. Data can be obtained from the E3SM website (https://e3sm.org/).
5. Download the default E3SM MOSART input data, including the river network data. Run the MOSART model to generate the river flow data. This data will be used as the boundary condition for the E3SM ELM model. An example data stream file is provided in the `data/` directory.
6. Run the PyE3SM to setup the E3SM ELM model with the hillslope capability, an example of the case setup is provided as `codes/k34/cases/h2sc/create_customized_h2sc_elm_case_2.py`. The hillslope data generated from the HexWatershed model is used to define the subgrid hillslope in the E3SM ELM model.
7. Run the E3SM ELM model. A default ELM spin up simulation is also recommended before running the hillslope simulation. See the above case setup file for more details.

## Reproduce my figures

Use the scripts found in the `codes/k34/analysis/plot/` directory to reproduce the figures used in this publication.