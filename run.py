# Input script for Mobility algorithm based on Lodhia & Clark (2022)
# Example for hydrogen
# Bhavik Harish Lodhia March 2024

import mobility

# See Appendix A in Hantschel & Kauerauf (2009) for description of lithological parameters for the following rock types

clastics = ["Sandstone","Sandstone-clayrich","Sandstone-claypoor", "Quartzite","Quartzite-quartz", \
            "Subarkose","Subarkose-quartz","Subarkose-clayrich","Subarkose-claypoor", "Subarkose-dolomite", \
            "Arkose","Arkose-quartzrich","Arkose-quartzpoor","Arkose-clayrich","Arkose-claypoor", "Arkose-dolomite", \
            "Siltstone","Siltstone-orgrich","Siltstone-TOC10","Siltstone-TOC2-3", "Conglomerate","Conglomerate-quartzite"]

carbonates = ["Limestone-OG","Limestone-WM","Micrite","Limestone-shaley","Limestone-orgrich","Limestone-TOC1-2", "Marl"]
dolomites = ["Dolomite","Dolomite-sandy","Dolomite-silty","Dolomite-org"]

# INPUTS
depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, 1.7, 1.8,1.9,2.0] # Values must be a list, i.e. [depth1, depth2, ...]
tsurf = 20. # Surface temperature. Value must be a float!
eos = "PR78" # Equation of state. See documentation for further information.
rock = "Sandstone" # Rock type. See documentation for list of rock types


mobility.Run.run(rock, depths, tsurf, eos, plot="on")