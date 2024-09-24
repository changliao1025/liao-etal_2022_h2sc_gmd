import os
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms

class RotationAwareAnnotation(mtext.Annotation):
    def __init__(self, s, xy, p, pa=None, ax=None, **kwargs):
        self.ax = ax or plt.gca()
        self.p = p
        if not pa:
            self.pa = xy
        self.calc_angle_data()
        kwargs.update(rotation_mode=kwargs.get("rotation_mode", "anchor"))
        mtext.Annotation.__init__(self, s, xy, **kwargs)
        self.set_transform(mtransforms.IdentityTransform())
        if 'clip_on' in kwargs:
            self.set_clip_path(self.ax.patch)
        self.ax._add_text(self)

    def calc_angle_data(self):
        ang = np.arctan2(self.p[1]-self.pa[1], self.p[0]-self.pa[0])
        self.angle_data = np.rad2deg(ang)

    def _get_rotation(self):
        return self.ax.transData.transform_angles(np.array((self.angle_data,)), 
                            np.array([self.pa[0], self.pa[1]]).reshape((1, 2)))[0]

    def _set_rotation(self, rotation):
        pass

    _rotation = property(_get_rotation, _set_rotation)

ninterval = 11
dSlope_surface = 40 / 850 
dArea = 3078 * 1.0e6 #m2
dHillslope_width = np.sqrt(dArea)
dLength_hillslope = dHillslope_width / 2
dLength_hillslope = 1367

dThickness_critical_zone = 80.0
dRatio0 = 0.5  #for bedrock slope

dRatio2 = 1.1 #above seepage
dRatio3 = 0.9 #below transition
nScenario = 19
#we assume that the soil thickness decreases at the elevation peak
dThickness_critical_zone_min = dThickness_critical_zone * dRatio0
dThickness_critical_zone_max = dThickness_critical_zone
#we assume the bedrock is almost paralle with surface
aDrop     = np.arange(5,81,5)
nScenario = len(aDrop)
aIncrease = np.arange(5, nScenario*5+5,5) #np.array([ 10.0, 30.0, 70.0, 140.0, 250.0])
nCell = 1

zwt = 44
topo= 61
for i in range(nCell):  
    #setup an elevation profile
    aElevation_profile = np.array([  38. ,  90,  95. ,  99,
                                    103. ,  107,  111. ,  116,
                                        120. ,  125,  418])
    #get min and max elevation
    dElevation_min = np.min(aElevation_profile)
    dElevation_max = np.max(aElevation_profile)

    dSlope_surface = 0.04 #40.0/850.0 #from in situ measurement
    
    #do we need to adjust elevation and length?
    #it is recommended to change the max 
    dElevation_max = dElevation_min + dLength_hillslope * dSlope_surface
    if (dElevation_max > 8848 ) :
      dElevation_max = 8848
      dSlope_surface = (dElevation_max - dElevation_min) / dLength_hillslope
      aElevation_profile[10] = dElevation_max
   
    
    
    #surface slope
    dElevation_difference_surface = dElevation_max - dElevation_min
    dSlope_surface = dElevation_difference_surface / dLength_hillslope
    #approximate the bedrock slope
    dElevation_difference_bedrock = ( dElevation_max - dThickness_critical_zone_min)  - ( dElevation_min - dThickness_critical_zone_max )
    dSlope_bedrock = dElevation_difference_bedrock / dLength_hillslope


    #set reference first   
    A1 = dSlope_surface
    dSlope_surface_radian = math.atan(A1)
    #now A1 line is known
    B1 = dElevation_min
    #reference transition water table slope
    dRatio1 = 1.0 - np.power(A1, 1.0)      
    dRatio1 = 0.25
    dSlope_water_table_reference= dSlope_surface * dRatio1
    A3 = dSlope_water_table_reference
    #now A3 line are known
    B3 = dElevation_min
    C3 = dLength_hillslope
    D3 = A3 * C3 + B3
    #set up range for different water table conditions
    dRange_left_without_seepage = dThickness_critical_zone
    dRange_right_without_seepage = dThickness_critical_zone + A3 * C3  
    #bedrock 
    A4 = dSlope_bedrock
    B4 = dElevation_min - dThickness_critical_zone
    
    
    #start drawing
    # start from here , we have another loop for WT dynamics
    fig = plt.figure()
    fig.set_figwidth( 20 )
    fig.set_figheight( 8 )
    plt.cla()
    ax = plt.axes()    
   
    #draw land surface
    x=np.array([0, dLength_hillslope])    
    y = x * A1 + B1
    ax.set_xlabel('Distance (m)', fontsize=10)
    ax.set_ylabel('Elevation (m)', fontsize=10)    
    fmt = '%.1E' 
    xticks = mtick.FormatStrFormatter(fmt)
    ax.xaxis.set_major_formatter(xticks)    
    ax.plot(x, y,'red',label='Surface')
    ax.set_xlim(x[0], x[1])

    ra = RotationAwareAnnotation("Downslope face", xy=(x[0],dElevation_min - dThickness_critical_zone), p=(x[0],dElevation_min), ax=ax,
                             xytext=(-10,0), textcoords="offset points", va="top")
    ra = RotationAwareAnnotation("Seepage face", xy=(x[0],dElevation_min), p=(x[1],y[1]), ax=ax,
                             xytext=(20,30), textcoords="offset points", va="top")

    #draw bedrock
    y = x * A4 + B4
    ax.plot(x,y, 'yellow',label='Bed rock')

    #a1_degree = np.arctan(A1) 
    #a2_degree = np.rad2deg( a1_degree)    
    #dummy0 = np.array((a2_degree,)) 
    #dummy1 = np.array( [ 0, dElevation_min ] ) 
    #dummy2 = dummy1.reshape((1, 2))
    #trans_angle = ax.transData.transform_angles(dummy0, dummy2, False )[0]
   
    #draw label 
    
    #ax.text(-0.03 * dLength_hillslope , dElevation_min - dThickness_critical_zone, "Downslope end", color ='blue', rotation = 90)
   
    #ax.text(0.06* dLength_hillslope, dElevation_min, "Seepage face", color ='green', rotation = trans_angle)
    ra = RotationAwareAnnotation("Bedrock", xy=(x[0],y[0]), p=(x[1],y[1]), ax=ax,
                             xytext=(2,-1), textcoords="offset points", va="top")

    iFlag_once = 1 

  
    #right range
    dRange_left_with_seepage = dElevation_difference_surface # dElevation_difference_bedrock
    dRange_right_with_seepage = dElevation_max - D3    

    #draw transition 
    G34 = (B3-B4) / (A4-A3)
    H34 = A4 * G34 + B4
    if(G34 < dLength_hillslope):
        x=np.array([0, G34])
        y = x * A3 + B3
        ax.plot(x,y, 'brown',label ='Water table transition')
        x= np.array([G34, dLength_hillslope])
        y = x * A3 + B3
        ax.plot(x,y, 'brown', linestyle='--' )
    else:
        x=np.array([0, dLength_hillslope])
        y = x * A3 + B3
        ax.plot(x,y, 'brown',label ='Water table transition')
        

    for j in range(nScenario):
        #water table below minimal elevation       

        
        dHeight1 = aDrop[j]       
        dElevation_water_table = dElevation_min - aDrop[j]
        B5 = dElevation_water_table    
        dDummy1 = dHeight1 / dThickness_critical_zone
        dDummy11 = np.power( dDummy1, dRatio3)
        if dDummy11 > 1.0:
            dDummy11 = 1.0
        else:
            print(dDummy1, dDummy11)

        dDummy2 = dRange_right_without_seepage * dDummy11
        D5 = D3 -  dDummy2
        F5 = dElevation_water_table
        dDummy3 = (D5 - F5) / dLength_hillslope      
        A5 = dDummy3   
        if A5 < 0:
            print("A5", A5)
        dSlope_watertable1 = math.atan( A5 )
        #with the slope known, we can also check the intersect length because it might intersect with the bed rock
        #intersect watertable below minimal elevation (A4 with A5)
        #y2 = A4 * x + B4 
        #y3 = A5 * x + B5         
        B5 = dElevation_min-dHeight1
        #x = (B5-B4) / (A4-A5)
        G45 = (B5-B4) / (A4-A5)
        H45 = A4 * G45 + B4
        dLength_watertable1 = G45

        #=========================================
        #water table above minimal elevation/seepage          
        #=========================================
        dHeight2 = aIncrease[j]
        if dHeight2 > dRange_left_with_seepage:
            dHeight2 = dRange_left_with_seepage
        dElevation_water_table = dElevation_min + dHeight2
        dDummy4 = dHeight2 / dRange_left_with_seepage        
        dDummy5 = np.power( dDummy4, dRatio2)        
        dDummy6 = dRange_right_with_seepage * dDummy5        
        D2 = D3 + dDummy6    
        H12 = dElevation_water_table        
        G12 = (H12 - B1)  /  A1
        C2 = dLength_hillslope
        if D2<H12:
            D2=H12

        A2 = ( D2 - H12)/( C2 - G12 )
        dSlope_watertable2 = math.atan( A2  )
        #=========================================
        #intersect above elevation (A4 A2)
        #=========================================
        #y2 = A4 * x + B4
        #y4 = A2 * x + B2
        #y41 = A2 * x41  + B2, where x41 and y41 are intersect between A1 and A2
        #B2 unknown, but it could be calculated from:        
        B2 = H12 - A2 * G12        
        G24 = (B2-B4) / (A4-A2)
        H24 = A4 * G24 + B4
        dLength_watertable2 = G24        
        #intersect between two watertable
        G25 = (B2-B5) / (A5-A2)
        H25 = A5 * G25 + B5  
        
        
        ci = 'C' + str(j)
        if iFlag_once == 1:
            x=np.array([0, G45])
            y = x * A5 + B5
            ax.plot(x,y, 'blue',label ='Water table without seepage')
            x= np.array([G45, dLength_hillslope])
            y = x * A5 + B5
            ax.plot(x,y, 'blue', linestyle='--' )
            
        else:
            x=np.array([0, G45])
            y = x * A5 + B5
            ax.plot(x,y, 'blue')
            x= np.array([G45, dLength_hillslope])
            y = x * A5 + B5
            ax.plot(x,y, 'blue', linestyle='--')
            
        #plot part
        x_1 = [0, G12]
        y_1 = [H12, H12]
        #intersect as well
        G24 = (B2-B4) / (A4-A2)
        H24 = A4 * G24 + B4
        if(G24 < dLength_hillslope):
            x_2 = np.array([G12, G24])
            y_2 = x_2 * A2 + B2
            x_22 = np.array([G24, dLength_hillslope])
            y_22 = x_22 * A2 + B2 
            if iFlag_once == 1:
                ax.plot(x_1,y_1, 'green', linestyle='--')
                ax.plot(x_2,y_2, 'green', label ='Water table with seepage')
                ax.plot(x_22,y_22, 'green', linestyle='--')
            
            else:
                ax.plot(x_1,y_1, 'green', linestyle='--')
                ax.plot(x_2,y_2, 'green')
                ax.plot(x_22,y_22, 'green', linestyle='--')
                
        else: 
            x_2 = np.array([G12, dLength_hillslope])
            y_2 = x_2 * A2 + B2
            if iFlag_once == 1:
                ax.plot(x_1,y_1, 'green', linestyle='--')
                ax.plot(x_2,y_2, 'green', label ='Water table with seepage')
            else:
                ax.plot(x_1,y_1, 'green', linestyle='--')
                ax.plot(x_2,y_2, 'green')
             

        
        iFlag_once = 0    
        
    
    ax.legend()
    print("=============")            

    sFilename = os.path.dirname(__file__) + '/slope_without_seepage' + '.png'
    plt.savefig(sFilename )

