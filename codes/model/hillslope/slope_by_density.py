
import numpy as np
gxr = 0.0003

dArea = 3078*1.0E6

dArea_half = dArea / 2.0

dChannel_length = dArea * gxr
hlen =  dArea / dChannel_length 

dmin=68
dmax=418

dslp = (dmax-dmin)/hlen

print ('H2SC hillslope issue', dslp)


