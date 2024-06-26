# hydrogen_mobility
Dr Bhavik Harish Lodhia
March 2024

Fluid mobility calculation using the Darcy Flow equation.
This module calculates density, buoyancy, viscosity and maximum vertical velocity
for hydrogen in various rock types using the method of Lodhia & Clark (2022).
Currently, only pure fluids are supported. The addition of other, mixed or multi-component
fluids is currently under development.

To calculate hydrogen fluid properties, please use the command: python run.py
Please see the detailed comments in run.py for information on how hydrogen fluid properties for different rock types and conditions may be calculated and plotting options.

Required modules:
pip install chemicals thermo pandas uncertainties matplotlib

Please cite use of this software as follows:
Lodhia, B.H. (2024) hydrogen_mobility, https://github.com/lodhia92/hydrogen_mobility, doi:10.5281/zenodo.10990921.

References:

Lodhia, B.H., Peeters, L. and Frery, E. (2024). A Review of the Migration of Hydrogen from the Planetary to Basin Scale. EarthArxiv, doi:10.31223/X5C11J.

Lodhia, B.H., Peeters, L. (2024). The migration of hydrogen in sedimentary basins, Australian Energy Producers Journal, doi:10.1071/EP23176

Lodhia, B.H., Clark, S.R. (2022). Computation of vertical fluid mobility of CO2, methane, hydrogen and hydrocarbons through sandstones and carbonates. Sci Rep 12, 10216. https://doi.org/10.1038/s41598-022-14234-6
