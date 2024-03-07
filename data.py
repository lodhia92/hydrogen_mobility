# Class to set up Mobility module data parameters
# Bhavik Harish Lodhia UNSW Sydney
# 17/05/2022

import pandas as pd  
from thermo import ChemicalConstantsPackage
from thermo.eos_mix import PR78MIX, TWUPRMIX, SRKMIX, TWUSRKMIX, APISRKMIX, RKMIX
from thermo.eos_mix import VDWMIX
#from chemicals.viscosity import Lorentz_Bray_Clarke, mu_IAPWS
#from chemicals.iapws import iapws95_rho, iapws97_rho
#from chemicals import Vm_to_rho
import numpy as np
  
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

    def PT(setting, depth, tsurf):
        # Assuming a geothermal gradient of 25 degrees/km, change as required
    
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
                                                 
###############################################################################
                                                 
class Fluid:

    def Name(fluid):    

        if fluid == "methane":
            name = ChemicalConstantsPackage(MWs=[16.043], names=['methane'],
                                               omegas=[0.0115], Pcs=[4.599e6], Tcs=[190.25])            
            vc = [0.0986/1000]
        elif fluid == "CO2":
            name = ChemicalConstantsPackage(MWs=[44.010], names=['CO2'], omegas=[0.2276], 
                                           Pcs=[7.382e6], Tcs=[304.19])
            vc = [0.0940/1000]
        elif fluid == "H2":
            name = ChemicalConstantsPackage(MWs=[2.016], names=['H2'], omegas=[-0.2150], 
                                          Pcs=[1.313e6], Tcs=[33.18])            
            vc = [0.0642/1000]
        elif fluid == "H2O":
            name = ChemicalConstantsPackage(MWs=[18.015], names=['H2O'], omegas=[0.3449], 
                                           Pcs=[22.055e6], Tcs=[647.13]) 
            vc = [0.0560/1000]
        elif fluid == "drygas":
            name = ChemicalConstantsPackage(MWs=[17.943], names=['dry_gas'], 
                                              omegas=[0.0221], Pcs=[4.850e6], Tcs=[197.42])            
            vc = [0.0977/1000]
        elif fluid == "wetgas":
            name = ChemicalConstantsPackage(MWs=[30.186], names=['wet_gas'], 
                                              omegas=[0.0624], Pcs=[4.801e6], Tcs=[272.40])            
            vc = [0.1345/1000] 
        elif fluid == "voil":
            name = ChemicalConstantsPackage(MWs=[53.135], names=['voloil'],
                                                   omegas=[0.1267],Pcs=[4.373e6], Tcs=[367.79])            
            vc = [0.1956/1000]
        elif fluid == "lightoil":
            name = ChemicalConstantsPackage(MWs=[48.439], names=['lightoil'],
                                                omegas=[0.1238],Pcs=[4.389e6], Tcs=[357.38])            
            vc = [0.1895/1000]
        elif fluid == "blackoil":
            name = ChemicalConstantsPackage(MWs=[90.072], names=['blackoil'],
                                                omegas=[0.2423],Pcs=[3.535e6], Tcs=[491.62])            
            vc = [0.3121/1000]
        elif fluid == "mediumoil":
            name = ChemicalConstantsPackage(MWs=[101.141], names=['mediumoil'],
                                                   omegas=[0.2606],Pcs=[3.459e6], Tcs=[511.92])     
            vc = [0.3302/1000]
        elif fluid == "heavyoil":
            name = ChemicalConstantsPackage(MWs=[127.492], names=['heavyoil'],
                                                omegas=[0.3290],Pcs=[3.041e6], Tcs=[568.42])       
            vc = [0.4080/1000]
        return name, vc









                                          