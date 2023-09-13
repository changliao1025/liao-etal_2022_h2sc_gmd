import numpy as np
nElevation_interval_in=11 
dElevation_water_table_in =  61.6-5
aElevation_profile_in = np.array([  38. ,  90,  95. ,  99,
                                    103. ,  107,  111. ,  116,
                                        120. ,  125,  418])
aSeepage_mask_out = np.zeros(nElevation_interval_in)
for iIndex in range(nElevation_interval_in-1):
    if aElevation_profile_in[iIndex] < dElevation_water_table_in:
        aSeepage_mask_out[iIndex] = 1

iIndex = int(np.sum(aSeepage_mask_out))
dLength_hillslope_in =27000

dummy_start = aElevation_profile_in[iIndex-1]
dummy_end = aElevation_profile_in[iIndex ]
dummy1 = (iIndex-1)/ nElevation_interval_in * dLength_hillslope_in
dummy2 = 1.0/nElevation_interval_in  * ( dElevation_water_table_in -dummy_start)/(dummy_end-dummy_start)* dLength_hillslope_in
dLength_seepage_out = dummy1  + dummy2

print(dummy1, dummy2, dLength_seepage_out)
    

       
