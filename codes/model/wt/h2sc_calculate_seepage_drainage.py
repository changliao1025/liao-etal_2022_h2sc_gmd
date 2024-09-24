import numpy as np

dSlope_surface_in = 40.0/850.0 #5.699999999999999E-003
dSlope_surface_in = 0.014
dSaturate_hydraulic_conductivity_in= 0.107872755330942

dArea = 3087.741 * 1000000
dLength_hillslope = np.sqrt(dArea) / 2.0* 1000 
#dLength_hillslope = 27000  * 1000 

dWidth_hillslope_in= dLength_hillslope * 2
dLength_seepage = 10 * 1000

dArea_hillslope_in = dWidth_hillslope_in*dLength_hillslope


#head pressure delta H, average  
dDummy1 = dLength_hillslope * dSlope_surface_in
#cross area A
dDummy2 = dLength_seepage  * dWidth_hillslope_in  
#q = k   A   delta H/ L

#delta h/L
dDummy3 = dSlope_surface_in

dDummy4 =  dSaturate_hydraulic_conductivity_in * dDummy2 * dDummy3 


dDummy5 = dArea_hillslope_in


dFlow_seepage_out = dSaturate_hydraulic_conductivity_in \
         * dLength_seepage * np.power(dSlope_surface_in,2) * dWidth_hillslope_in / dArea_hillslope_in


print(dFlow_seepage_out)