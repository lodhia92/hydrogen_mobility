# Class to set up Mobility module data parameters
# Bhavik Harish Lodhia

import pandas as pd  
from thermo import ChemicalConstantsPackage
  
###############################################################################   


class Data:
    
    def Params():

        header_list = ["Rock","ak","phi0","phi1","phi2","k0","k1","k2"]     
        df = pd.read_csv("data/permeability.csv", names=header_list)        
#        x = df[['Rock', 'ak', 'phi0','k0','phi1','k1','phi2','k2']]
        compaction_headers = ["Rock","phi0","athy_k_km","athy_k_MPa","Cmax","Cmin", \
                              "ka","kb","phi"]
        compaction = pd.read_csv("data/compaction.csv", names=compaction_headers)        
        athy = compaction[['Rock', 'phi0', 'athy_k_km']]   
        
        return df, athy, compaction

    def Name(rock):

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

    def PT(depth, tsurf):
        # Assuming a geothermal gradient of 25 degrees/km, change as required
    
        temp = (depth*25) + tsurf  # depth in km for consistency with other functions
                                   # temperature in degrees celsius, convert to K
                                   # temperature at surface set to 20 degrees

        m = 0.5 # Typical PT gradient for normal geological conditions = 0.5 MPa/k

    
        # calculate graph intercept from tsurf
        C = -1*tsurf*m    
    
        pressure = m*temp + C
   
        # pressure !<0, set surface pressure = 1 atm = 0.10325 MPa     
        if pressure < 0.101325:
            pressure = 0.101325
        
        
        return temp + 273.15, pressure  # return temperature in kelvin, 
                                                 # pressure in MPa   
                                                 
###############################################################################

class Fluid:
                                                     
    def Name(fluid):    

       
        if fluid == "H2":
            name = ChemicalConstantsPackage(MWs=[2.016], names=['H2'], omegas=[-0.2150], 
                                          Pcs=[1.313e6], Tcs=[33.18])            
            vc = [0.0642/1000]
        elif fluid == "H2O":
            name = ChemicalConstantsPackage(MWs=[18.015], names=['H2O'], omegas=[0.3449], 
                                           Pcs=[22.055e6], Tcs=[647.13]) 
        
        return name, vc









                                          