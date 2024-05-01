# Input script for Mobility algorithm based on Lodhia & Clark (2022)
# Example for hydrogen
# Bhavik Harish Lodhia 06/03/2024

import mobility

sandstone = ["Sandstone","Sandstone-clayrich","Sandstone-claypoor"]
quartzite = ["Quartzite","Quartzite-quartz"]
arkose = ["Subarkose","Subarkose-quartz","Subarkose-clayrich","Subarkose-claypoor", "Subarkose-dolomite", \
            "Arkose","Arkose-quartzrich","Arkose-quartzpoor","Arkose-clayrich","Arkose-claypoor", "Arkose-dolomite"]
siltstone = ["Siltstone","Siltstone-orgrich","Siltstone-TOC10","Siltstone-TOC2-3"]
conglomerate = ["Conglomerate","Conglomerate-quartzite"]

limestone = ["Limestone-OG","Limestone-WM","Micrite","Limestone-shaley","Limestone-orgrich","Limestone-TOC1-2"]
carbonates = ["Marl","Dolomite","Dolomite-sandy","Dolomite-silty","Dolomite-org"]

# INPUTS

variable = "vmax" # Choose from vmax, mobility, buoyancy, density, viscosity
setting = 3 # Geological settings: 1 = 1 = cool overpressured, 2.5 MPa/K, 2 = overpressured, 1.0 MPa/K, 3 = typical, 0.5 MPa/k
            # 4 = hydrostatic, 0.3 MPa/K, 5 = hot hydrostatic, 0.1 MPa/K (see Lodhia & Clark 2022 for more information)
            # Default choice for setting is normal geological conditions.
depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, \
          1.7, 1.8,1.9,2.0] # Values must be a list, i.e. [depth1, depth2, ...]
tsurf = 20. # Surface temperature. Value must be a float!
eos = "PR78" # Equation of state. See documentation for further information.
rock = "Sandstone" # Rock type. See documentation for list of rock types
fluid = "H2" # Example for Hydrogen

mobility.Run.run(fluid,rock,depths,tsurf,setting,eos,output="on",plot=variable,save="false") # Requires fluid, rock, depths (must be a list), surface temperature,
                                                                                                # EOS, output ("on" or "off") displayed in terminal, plot ("name of variable")
                                                                                                # save ("on" will save output as csv file and plot as png file, anything
                                                                                                # else will not save)