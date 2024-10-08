

# Simulation configurations

# K34 cases

| Case |Model| Length | Parameters | Note |  Purpose |
|---|-|--|---------|-----------------|--|
| 0 | Default WT-based | 90 years (1890-1979) |Default | A simple custmoized ELM simulation | For ELM spinup |


| Case |Model| Length | Parameters | Note |  Purpose |
|---|-|--|---------|-----------------|--|
| 1 | Default WT-based |30 year (1980-2009)| Default | A custmoized ELM simulation, use Case 0 as initial | For default transient |
| 2 | HLG-based |30 year (1980-2009) |Default (slope + anisotropy ratio) | A custmoized ELM-DROF simulation, use Case 0 as initial  | Baseline hillslope-based model |
| 3 | HLG-based |30 year (1980-2009) |Actual hill slope | A custmoized ELM-DROF simulation, use Case 0 as initial  | The role of slope in hillslope-based model |
| 4 | HLG-based |30 year (1980-2009) |Use surface slope as WT slope | A custmoized ELM-DROF simulation, use Case 0 as initial  | The role of slope in hillslope-based model |
| 5 | HLG-based |30 year (1980-2009) |Anisotropy ratio as 1.0 | A custmoized ELM-DROF simulation, use Case 0 as initial  | The role of anisotropy ratio in hillslope-based model |
| 6 | HLG-based |30 year (1980-2009) |Gage height time invariant 0.5 m | A custmoized ELM-DROF simulation, use Case 0 as initial  | The role of gage height in hillslope-based model |
| 7 | HLG-based |30 year (1980-2009) |Gage height time invariant 10.0 m | A custmoized ELM-DROF simulation, use Case 0 as initial  | The role of gage height in hillslope-based model |


sDate = 20230401 


## Analysis

|  | Cases used |   Purpose |
|---|-----------|-----------------|
| 1 | 1,2,3 |  The role of slope in subsurface drainage  |
| 2 | 1,2,3  |  The role of slope in WTD  |
| 3 | 1,2,4  |  The role of anisotropy ratio in subsurface drainage  |
| 4 | 1,2,4  |  The role of anisotropy ratio in WTD  |
| 5 | 2,5,6 |  The role of gage height in subsurface drainage  |
| 6 | 2,5,6 |  The role of gage height in WTD  |
| 7 | 2,5,6 |  The role of gage height in WTD  |
| 8 | 2,5,6 |  The role of surface slope in WTD  |