#H2SC hillslope issue 1 511.961663509409 240000.000000000
#H2SC hillslope issue 1 511.961663509409 7838.07371878634
#H2SC hillslope issue 0.304000000000000 122870799.242258
#240000.000000000 1250.00000000000 245741598.484516
import numpy as np

dArea = 3078*1.0E6
dArea_half = dArea / 2.0
hlen = np.sqrt(dArea) / 2.0# dArea / dChannel_length 

dmin=68
dmax=418

dslp = (dmax-dmin)/hlen


print ('H2SC hillslope issue', dslp)


