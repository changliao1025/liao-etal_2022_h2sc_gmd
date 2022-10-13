#https://books.gw-project.org/hydrogeologic-properties-of-earth-materials-and-principles-of-groundwater-flow/chapter/equation-derivation-for-equivalent-k-and-a-4-layer-application/
import numpy as np

aLayer_thickness = np.arange(1, 6)

aHK_layer_vertical = np.arange(1, 6) *0.1

aAnisotropy = np.arange(1, 6) * 2

d=np.sum(aLayer_thickness)

aHK_layer_horizontal = aHK_layer_vertical * aAnisotropy

#method 1

Kx= np.sum( aHK_layer_horizontal* aLayer_thickness)/d

print(Kx)

#method 2

Kz = d / np.sum( aLayer_thickness/ aHK_layer_horizontal)

print(Kz)