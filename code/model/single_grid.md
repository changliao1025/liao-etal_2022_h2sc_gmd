To set up a single grid case, we need to modify the workflow in the following steps:

1. Generate a time series MOSART dataset, which includes water level. This will be used as the river conditions. 
2. Add the river condition as a forcing data in ELM
3. Modify ELM to read the river condition
4. Prepare a surface dataset that has all the hillslope information

Steps:

1. Change compset setting:
https://github.com/E3SM-Project/E3SM/compare/master...donghuix:E3SM:donghuix/ocn-lnd/one-way-coupling#diff-114250fbefb5fc4258e488a61018e7a94e9e9693374e9221b6a928977a4c9a42


https://github.com/donghuix/E3SM/tree/donghuix/ocn-lnd/one-way-coupling