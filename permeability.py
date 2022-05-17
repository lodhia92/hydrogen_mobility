# Module to calculate rock (intrinsic) permeability using measured points method
# Bhavik Harish Lodhia UNSW Sydney
# 28/06/21

import pandas as pd
import math
from uncertainties import ufloat

###############################################################################   

# Assign column names
header_list = ["Rock","ak","phi0","phi1","phi2","k0","k1","k2"] 

# Load .csv file with columns
df = pd.read_csv("permeability.csv", names=header_list)        

# Reorder and create dataframe from .csv file using Pandas
x = df[['Rock', 'ak', 'phi0','k0','phi1','k1','phi2','k2']]
# print(x)

compaction_headers = ["Rock","phi0","athy_k_km","athy_k_MPa","Cmax","Cmin", \
                      "ka","kb","phi"] # Assign column names

# Load .csv file with columns
compaction = pd.read_csv("compaction.csv", names=compaction_headers)        

# Reorder and create dataframe from .csv file using Pandas
athy = compaction[['Rock', 'phi0', 'athy_k_km']]    

# print(athy)

###############################################################################   
def name(rock):

    # Assign row number for each lithology
    if rock == "Limestone-OG":
        row = 0
    elif rock == "Limestone-WM":
        row = 1
    elif rock == "Micrite":
        row = 2
    elif rock == "Limestone-shaley":
        row = 3   
    elif rock == "Limestone-orgrich":
        row = 4   
    elif rock == "Limestone-TOC1-2":
        row = 5        
    elif rock == "Limestone-TOC10":
        row = 6               
    elif rock == "Marl":
        row = 7              
    elif rock == "Dolomite":
        row = 8          
    elif rock == "Dolomite-sandy":
        row = 9         
    elif rock == "Dolomite-silty":
        row = 10           
    elif rock == "Dolomite-org":
        row = 11
    elif rock == "Chalk":
        row = 12   
    elif rock == "Chalk-calcite95":
        row = 13           
    elif rock == "Chalk-calcite75":
        row = 14            
    elif rock == "Chalk-calcite40":
        row = 15   
    elif rock == "Coal":
        row = 16 
    elif rock == "Coal-impure":
        row = 17 
    elif rock == "Coal-silty":
        row = 18            
    elif rock == "Sandstone":
        row = 19
    elif rock == "Sandstone-clayrich":
        row = 20
    elif rock == "Sandstone-claypoor":
        row = 21
    elif rock == "Quartzite":
        row = 22
    elif rock == "Quartzite_quartz":
        row = 23
    elif rock == "Subarkose":
        row = 24
    elif rock == "Subarkose-quartz":
        row = 25
    elif rock == "Subarkose-clayrich":
        row = 26
    elif rock == "Subarkose-claypoor":
        row = 27
    elif rock == "Subarkose-dolomite":
        row = 28
    elif rock == "Arkose":
        row = 29
    elif rock == "Arkose-quartzrich":
        row = 30
    elif rock == "Arkose-quartzpoor":
        row = 31
    elif rock == "Arkose-clayrich":
        row = 32
    elif rock == "Arkose-claypoor":
        row = 33
    elif rock == "Arkose-dolomite":
        row = 34
    elif rock == "Wacke":
        row = 35     
    elif rock == "Shale":
        row = 36         
    elif rock == "Shale-orglean":
        row = 37          
    elif rock == "Shale-sandy":
        row = 38          
    elif rock == "Shale-silty":
        row = 39          
    elif rock == "Shale-silicious":
        row = 40        
    elif rock == "Shale-opalCT":
        row = 41  
    elif rock == "Shale-black":
        row = 42  
    elif rock == "Shale-orgrich":
        row = 43  
    elif rock == "Shale-TOC3":
        row = 44  
    elif rock == "Shale-TOC8":
        row = 45  
    elif rock == "Shale-TOC20":
        row = 46  
    elif rock == "Siltstone":
        row = 47
    elif rock == "Siltstone-orgrich":
        row = 48
    elif rock == "Siltstone-TOC10":
        row = 49
    elif rock == "Siltstone-TOC2-3":
        row = 50
    elif rock == "Conglomerate":
        row = 51
    elif rock == "Conglomerate-quartzite":
        row = 52
    elif rock == "Tuff-felsic":
        row = 53        
    elif rock == "Tuff-basaltic":
        row = 54  

    return row          
 
###############################################################################   
    
def k(method, rock, porosity):     
    
    row = name(rock)    
    
    # Assign variables from dataframe depending on lithology
    ak = df.values[row,1]
    phi0 = df.values[row,2]
    phi1 = df.values[row,3]
    phi2 = df.values[row,4]
    k0 = df.values[row,5]
    k1 = df.values[row,6]
    k2 = df.values[row,7]
    
    # Calculate permeability from linear multipoint function
    if method == "multipoint" and porosity < phi1:
        x = abs(k1 - k0)/(phi1 - phi0)*porosity + k0
        
    elif method == "multipoint" and porosity >= phi1 and porosity <= phi2:
        x = (k2 - k1)/(phi2 - phi1)*porosity + \
            (k2 - (k2 - k1)*(phi2 - phi0)/(phi2 - phi1))
        
    elif method == "multipoint" and porosity > phi2:
        x = k2
    
    elif rock == "Sandstone" or rock == "Sandstone-clayrich" or  \
        rock == "Sandstone-claypoor" and method == "KC":
        
        S0 =  10**6
        phi = porosity - (3.1*10**-10)*S0
        if phi <= 0.1:
            y = (2*10**16)*10*(phi/(S0**2 * (1-phi)**2))
        elif phi > 0.1:
            y = (2*10**14)*10*(phi**3/(S0**2 * (1-phi)**2))  
        
    elif rock == "Siltstone" or rock == "Siltstone-orgrich" or  \
        rock == "Siltstone-TOC10" or rock == "Siltstone-TOC2-3" and method == "KC":
            
        S0 = 10**7
        phi = porosity - (3.1*10**-10)*S0 
        if phi <= 0.1:            
            y = (2*10**16)*0.5*(phi**5/(S0**2 * (1-phi)**2))
        elif phi > 0.1:
            y = (2*10**14)*0.5*(phi**3/(S0**2 * (1-phi)**2))                
 
    elif rock == "Shale" or rock == "Shale-sandy" or rock == "Shale-silicious"  \
        or rock == "Shale-silty" or rock == "Shale-opalCT" or rock == "Shale-black" \
        or rock == "Shale-orgrich" or rock == "Shale-TOC3" or rock == "Shale-TOC8" \
        or rock == "Shale-TOC20" and method == "KC":
        S0 = 10**8
        phi = porosity - (3.1*10**-10)*S0 
        if phi <= 0.1:            
            y = (2*10**16)*0.01*(phi**5/(S0**2 * (1-phi)**2))
        elif phi > 0.1:
            y = (2*10**14)*0.01*(phi**3/(S0**2 * (1-phi)**2))


    if method == "multipoint":

        khv = 10**x # Calculate permeability in mD
        kv = khv*1 # vertical permeability = khv * upscaling
        kh = ak*khv*50 # horizonal permeability = vertical permiability * ak * upscaling
    
        kv = math.log(kv, 10) # math.log returns log_e unless specified
        kh = math.log(kh, 10)
        # print(porosity, kv,kh)
        
    elif method == "KC":
        khv = y
        kv = khv*1
        kh = ak*khv
        kv = math.log(kv, 10) # math.log returns log_e unless specified
        kh = math.log(kh, 10)
        # print(porosity, kv,kh)
    
    return kv, kh # vertical and horizontal permeability in log[mD]    
    
    # print(porosity, x)
    # return x


# for phi in range(0,101,1):
#     por = k("multipoint", "Micrite", phi)
#     # print(phi, por)


###############################################################################    

def krp(rock, equation, regime, S):
    
    # Relative permeability calculation using quadratic formula from Hantschel
    # (2009), Ringrose and Corbett (1994) methods
    
    if rock == "Sandstone" or "Siltstone" or "Conglomerate" or "Limestone" \
        or "Marl" or "Dolomite":
            Swc = 5/100
            Sgc = 0.00
            Soc = 0.10/100
            
    elif rock == "Shale":
            Swc = 5/100
            Sgc = 0.00
            Soc = 1.00/100     
            
    if equation == "Quadratic" and regime == "WL":        # water-liquid phase
        Swe = (S-Swc)/(1 - Swc - Soc)           # normalised water saturation
        krw = 0.4*(Swe)**2                      # water relative permeability
        krow = 1 - 1.8*(Swe) + 0.8*(Swe)**2     # oil-water relative permeability
        # print(str('Quadratic'), S, str('krw'), krw, str('krow'), krow)
        return S, krw, krow
    
    elif equation == "Quadratic"  and regime == "VL":    # vapour-liquid phase
        Sge = (S-Sgc)/(1 - Swc - Sgc) # normalised gas saturation for krg
        Sgoe = S/(1 - Swc)            # normalised gas-oil saturation for krog
        krg = 0.4*(Sge)**2            # relative gas permeability
        krog = 1 - 1.8*(Sgoe) + 0.8*(Sgoe)**2   # relative oil-gas permeability
        # print(str('Quadratic'), S, str('krg'), krg, str('krog'), krog)        
        return S, krg, krog

    elif equation == "Ringrose" and regime == "WL":
        Swe = (S-Swc)/(1 - Swc - Soc)           # normalised water saturation
        krw = 0.3*Swe**3
        krow = 0.85*(1 - Swe)**3
        # print(str('Ringrose'), S, str('krw'), krw, str('krow'), krow)        
        return S, krw, krow

    elif equation == "Ringrose"  and regime == "VL":
        print(str('Rongrose and Corbett 1994 method can only be used for water-liquid components'))
    
    
# for Sw in np.linspace(0, 1, 101):
#     krp("Shale", "Quadratic", "VL", Sw) 

###############################################################################   
    
def porosity(rock, depth):
    # Athy (1930) formula with depth in METRES and porosity given from 0 to 1
    
    row = name(rock)
    
    # Assign variables from dataframe depending on lithology
    dpor = athy.values[row,1]  # depositional porosity
    athyk = athy.values[row,2] # Athy compaction wavelength in km, *1000 to use metres
    
    porosity = (dpor*math.exp(-depth/athyk))/100
    return porosity

###############################################################################   

def Swi(rock, porosity):
    # Calculates connate water saturation from porosity using Holmes (2009) equation. 
    # Also calculates uncertainties using mid-point values for each constant
    # For sandstones, 0.02 < C < 0.1
    # For carbonates, 0.005 < C < 0.06
    # For all rock types, 0.8 < Q < 1.3
    
    # porosity**(Q) * Swi = constant
    
    Q = ufloat(1.05, 0.25) # 0.8 < Q < 1.3 for sandstones and cbates (Holmes 2009)

    
    if rock == "Sandstone":  
        dpor = 0.4101  # depositional porosity of typical sandstone (Hantschel 2009)
                       # add on 0.001 to dpor so dpor value is included in rounding
        C = ufloat(0.06, 0.04)
    elif rock == "Carbonate":
        dpor = 0.5101  # depositional porosity of micrite (Hantschel 2009)
                       # add on 0.001 to dpor so dpor value is included in rounding
        C = ufloat(0.0325,0.0275)        
    
    # Only return values for porosities less than depositional porosity!
    # Swi cannot be > 1, this is impossible. Errors increase for low porosities,
    # so I have set a condition that does not allow Swi > 1 +/- 1
    
    if porosity <= dpor:
        Sw = C/(porosity**Q)
        if Sw > 1:
            Sw = ufloat(1,1)
        return porosity, Sw
    elif porosity > dpor:
        pass

###############################################################################   

def SwiZ(rock, depth):
    
    # Calculates connate water saturation from porosity using Holmes (2009) 
    # equation and porosity from depth using Athy (1930)
    # Also calculates uncertainties using mid-point values for each constant
    # For sandstones, 0.02 < C < 0.1    - approximated for all clastics
    # For carbonates, 0.005 < C < 0.06  - approximated for all carbonates
    # For all rock types, 0.8 < Q < 1.3
    
    row = name(rock)
    
    # Assign variables from dataframe depending on lithology
    dpor = athy.values[row,1]  # depositional porosity
    athyk = athy.values[row,2] # Athy compaction wavelength, *1000 for metres
    
    porosity = (dpor*math.exp(-depth/athyk))/100
    # return porosity

    Q = ufloat(1.05, 0.25) # 0.8 < Q < 1.3 for sandstones and cbates (Holmes 2009)
    
    # Holmes (2009) parameters
    
        
    # for SANDSTONES  (I'm using this for all clastics here!)
    if rock == "Sandstone" or rock == "Sandstone-clayrich" \
    or rock == "Sandstone-claypoor" or rock == "Quartzite" \
    or rock == "Quartzite-quartz" or rock == "Subarkose" \
    or rock == "Subarkose-quartz" or rock == "Subarkose-clayrich" \
    or rock == "Subarkose-claypoor" or rock == "Subarkose-dolomite" \
    or rock == "Arkose" or rock == "Arkose-quartzrich" \
    or rock == "Arkose-quartzpoor" or rock == "Arkose-clayrich" \
    or rock == "Arkose-claypoor" or rock == "Arkose-dolomite" \
    or rock == "Siltstone" or rock == "Siltstone-orgrich" \
    or rock == "Siltstone-TOC10" or rock == "Siltstone-TOC2-3" \
    or rock == "Conglomerate" or rock == "Conglomerate-quartzite":

        C = ufloat(0.06, 0.04)
    
    
    # for CARBONATES
    elif rock == "Limestone-OG" or rock ==  "Limestone-WM" or rock ==  "Micrite" \
    or rock == "Limestone-shaley" or rock == "Limestone-orgrich" \
    or rock == "Limestone-TOC1-2" or rock == "Limestone-TOC10" \
    or rock == "Marl" or rock == "Dolomite" or rock == "Dolomite-sandy" \
    or rock == "Dolomite-silty" or rock == "Dolomite-org":
        
        # C = ufloat(0.0325,0.0275)
        C = ufloat(0.035,0.025)

    # Holmes (2009) equation
    Sw = C/(porosity**Q)
    if Sw > 1:
        Sw = ufloat(1,1)         # set Swi = 1.0 +/-1 for low porosities
    return depth, porosity, Sw
    # elif porosity > dpor:
    #     pass          
###############################################################################   

def PT(setting, depth, tsurf):
    # Assuming a geothermal gradient of 25-30 degrees/km
    
    temp = (depth*25) + tsurf  # depth in km for consistency with other functions
                               # temperature in degrees celsius, convert to K
                               # temperature at surface set to 20 degrees
    
    # PT gradients for each geological setting
    if setting == 1:         # 1 = cool overpressured, 2.5 MPa/K
        m = 2.5
    elif setting ==2:        # 2 = overpressured, 1.0 MPa/K
        m = 1.0
    elif setting == 3:       # 3 = typical, 0.5 MPa/k
        m = 0.5
    elif setting == 4:       # 4 = hydrostatic, 0.3 MPa/K
        m = 0.3
    elif setting == 5:       # 5 = hot hydrostatic, 0.1 MPa/K
        m = 0.1
    
    # calculate graph intercept from tsurf
    C = -1*tsurf*m    
    
    pressure = m*temp + C
   
    # pressure !<0, set surface pressure = 1 atm = 0.10325 MPa     
    if pressure < 0.101325:
        pressure = 0.101325
        
        
    return setting, temp + 273.15, pressure  # return temperature in kelvin, 
                                             # pressure in MPa
    

# tsurf = 40.
# for setting in range(1,6,1):
#     for z in range(0,10,1):
#         ptpath = PT(setting, z, tsurf)
#         print(ptpath)
        
    
    
# for z in np.linspace(0,6,61):
#     Sw = SwiZ("Micrite", z)
#     print(Sw)

# for por in np.linspace(0.01, 1.0, 100):
#     Sw = Swi("Carbonate", por)
#     print(Sw)   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    