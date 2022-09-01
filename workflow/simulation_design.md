

# Simulation configurations

## Cases 

| Case |Model| Length | Parameters | Note |  Purpose |
|---|-|--|---------|-----------------|--|
| 1 | Default WT-based | 90 years (1890-1979) |Default | A simple custmoized ELM simulation | For ELM spinup |
| 2 | Default WT-based |30 year (1980-2009)| Default | A custmoized ELM simulation, use Case 1 as initial | For default transient |
| 3 | H2SC-based |30 year (1980-2009) |Default (slope + anisotropy ratio) | A custmoized ELM-MOSART simulation, use Case 1 as initial  | Baseline hillslope-based model |
| 4 | H2SC-based |30 year (1980-2009) |WT slope as surface slope| A custmoized ELM-MOSART simulation, use Case 1 as initial  | The role of slope in hillslope-based model |
| 5 | H2SC-based |30 year (1980-2009) |Anisotropy ratio as 1.0 | A custmoized ELM-MOSART simulation, use Case 1 as initial  | The role of anisotropy ratio in hillslope-based model |
| 6 | H2SC-based |30 year (1980-2009) |Gage height time invariant 1.0 m | A custmoized ELM-MOSART simulation, use Case 1 as initial  | The role of gage height in hillslope-based model |
| 7 | H2SC-based |30 year (1980-2009) |Gage height time invariant 10.0 m | A custmoized ELM-MOSART simulation, use Case 1 as initial  | The role of gage height in hillslope-based model |




## Analysis

|  | Cases used |   Purpose |
|---|-----------|-----------------|
| 1 | 2, 3, 4 |  The role of slope in subsurface drainage  |
| 2 | 2, 3, 4 |  The role of slope in WTD  |
| 3 | 2, 3, 5 |  The role of anisotropy ratio in subsurface drainage  |
| 4 | 2, 3, 5 |  The role of anisotropy ratio in WTD  |
| 5 | 3, 6, 7 |  The role of gage height in subsurface drainage  |
| 6 | 3, 6, 7 |  The role of gage height in WTD  |
| 7 | 3, 6, 7 |  The role of gage height in WTD  |