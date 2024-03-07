# Input script for Mobility algorithm based on Lodhia & Clark (2022)
# Example for hydrogen, change as needed
# Bhavik Harish Lodhia 06/03/2024

import mobility
import data
import viscosity

depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, \
          1.7, 1.8,1.9,2.0] # Values must be a list, i.e. [depth1, depth2, ...]

#depth = [0.0] # Value must be a list, i.e. [depth]
tsurf = 20. # Surface temperature. Value must be a float!
eos = "PR78" # Equation of state. See documentation for further information.
rock = "Sandstone" # Rock type. See documentation for list of rock types
fluid = "methane" # Fluid can be H2, CO2, methane, dry_gas, wet_gas, lightoil, mediumoil, heavyoil, voil (volatile oil), blackoil

mobility.Run.run(fluid,rock,depths,tsurf,eos,output="on",plot="on") # Requires fluid, rock, depths (must be a list), surface temperature,
                                                                     # EOS, output ("on" or "off") displayed in terminal, plot ("on" or "off")