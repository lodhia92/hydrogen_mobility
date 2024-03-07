# Mobility
Dr Bhavik Harish Lodhia
March 2024

Fluid mobility calculation using the Darcy Flow equation.
This module calculates density, buoyancy, viscosity and maximum vertical velocity
for user-defined fluids in various rock types using the method of Lodhia & Clark (2022).
Currently, only pure fluids are supported. The addition of mixed or multi-component
fluids is currently under development.

Required modules:
pip install chemicals thermo pandas uncertainties matplotlib

References:

Lodhia, B.H., Peeters, L. (2024). The migration of hydrogen in sedimentary basins, Australian Energy Producers Journal, doi:10.1071/EP23176
Lodhia, B.H., Peeters, L., Frery, E. (2024). A Review of The Migration of Hydrogen from the Planetary to Basin Scale, UNDER REVIEW, Journal of Geophysical Research: Solid Earth.
Lodhia, B.H., Clark, S.R. (2022). Computation of vertical fluid mobility of CO2, methane, hydrogen and hydrocarbons through sandstones and carbonates. Sci Rep 12, 10216. https://doi.org/10.1038/s41598-022-14234-6
