# Module to calculate rock (intrinsic) permeability using measured points method
# Bhavik Harish Lodhia

import math
from uncertainties import ufloat
import data as data    
 
###############################################################################   

class Permeability:
    
    def k(rock, porosity):   
        
        df, athy, compaction = data.Data.Params()    
        row = data.Data.Name(rock)    
    
        # Assign variables from dataframe depending on lithology
        ak = df.values[row,1]
        phi0 = df.values[row,2]
        phi1 = df.values[row,3]
        phi2 = df.values[row,4]
        k0 = df.values[row,5]
        k1 = df.values[row,6]
        k2 = df.values[row,7]
    
        # Calculate permeability from linear multipoint function
        if porosity < phi1:
            x = abs(k1 - k0)/(phi1 - phi0)*porosity + k0
        
        elif porosity >= phi1 and porosity <= phi2:
            x = (k2 - k1)/(phi2 - phi1)*porosity + \
                (k2 - (k2 - k1)*(phi2 - phi0)/(phi2 - phi1))
        
        elif porosity > phi2:
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
 
        khv = 10**x # Calculate permeability in mD
        kv = khv*1 # vertical permeability = khv * upscaling
        kh = ak*khv*50 # horizonal permeability = vertical permiability * ak * upscaling
    
        kv = math.log(kv, 10) # math.log returns log_e unless specified
        kh = math.log(kh, 10)
        # print(porosity, kv,kh)
    
        return kv, kh # vertical and horizontal permeability in log[mD]    

###############################################################################    

    def krp(rock, regime, S):
    
        # Relative permeability calculation using quadratic formula from Hantschel
        # (2009)
    
        if rock == "Sandstone" or "Siltstone" or "Conglomerate" or "Limestone" \
            or "Marl" or "Dolomite":
            Swc = 5/100
            Sgc = 0.00
            Soc = 0.10/100
            
        elif rock == "Shale":
            Swc = 5/100
            Sgc = 0.00
            Soc = 1.00/100     
            
        if regime == "WL":        # water-liquid phase
            Swe = (S-Swc)/(1 - Swc - Soc)           # normalised water saturation
            krw = 0.4*(Swe)**2                      # water relative permeability
            krow = 1 - 1.8*(Swe) + 0.8*(Swe)**2     # oil-water relative permeability
            return S, krw, krow
    
        elif regime == "VL":    # vapour-liquid phase
            Sge = (S-Sgc)/(1 - Swc - Sgc) # normalised gas saturation for krg
            Sgoe = S/(1 - Swc)            # normalised gas-oil saturation for krog
            krg = 0.4*(Sge)**2            # relative gas permeability
            krog = 1 - 1.8*(Sgoe) + 0.8*(Sgoe)**2   # relative oil-gas permeability
     
            return S, krg, krog
    
###############################################################################   
    
    def Porosity(rock, depth):
        # Athy (1930) formula with depth in METRES and porosity given from 0 to 1
        athy = data.Data.Params()[1]
    
        row = data.Data.Name(rock)
    
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
    
        row = data.Data.Name(rock)
        athy = data.Data.Params()[1]
    
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
            
            C = ufloat(0.035,0.025)

        # Holmes (2009) equation
        Sw = C/(porosity**Q)
        if Sw > 1:
            Sw = ufloat(1,1)         # set Swi = 1.0 +/-1 for low porosities
        return depth, porosity, Sw
     
###############################################################################       