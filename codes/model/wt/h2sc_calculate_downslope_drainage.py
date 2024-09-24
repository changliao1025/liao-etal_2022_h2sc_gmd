import numpy as np
import math

dDummy_slope = 5.699999999999999E-003    #ratio. ,not degree or radian
dDummy_slope = 0.01348928 #ratio
sand = 10
xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) )
dSaturate_hydraulic_conductivity_in= 1 * 1.0E-3 #mm/s
dSaturate_hydraulic_conductivity_in = 0.107872755330942 #mm/s

dSaturate_hydraulic_conductivity_in = 0.1500848 #mm/s

dHeight_below_river_in=42*1000
dHeight_below_river_in=100*1000 #mm

dArea = 3078*1.0E6 #m2
dArea_half = dArea / 2.0  #m2
dResolution_grid_in = np.sqrt(dArea) #m
dLength_hillslope =  dResolution_grid_in / 2.0 #m
    
dDummy0 =  dLength_hillslope * 1000 #mm
#head pressure delta H, average  
dDummy1 = dDummy0 * (dDummy_slope)
#cross area A
dDummy2 = dHeight_below_river_in  * (dResolution_grid_in  * 1000) #m2
#q = k   A   delta H/ L

#delta h/L
dDummy3 = dDummy1 / dDummy0

dDummy4 =  dSaturate_hydraulic_conductivity_in * dDummy2 * dDummy3 #mm/s * mm2 * mm/mm=mm3/s


#drainage area, normalization
dDummy5 = dArea*0.5 * 1000000 #mm2

dFlow_downslope_out = dDummy4 /  dDummy5 #mm3/s / mm2 = m/s

dFlow_downslope_out = dFlow_downslope_out
#unit?

print(dFlow_downslope_out)


