# Calculate fluid velocities as a function of rock type, pressure and
# temperature
# By Bhavik Harish Lodhia

import data
import permeability
import viscosity
from uncertainties import ufloat
import pandas as pd
import matplotlib.pyplot as plt

class Mobility:
    
    # Calculate fluid mobility as a function of depth
    def Mobility(Fluid, rock, depth, tsurf, EOS):
    
        name, Vc = data.Fluid.Name(Fluid)
            
        # calculate porosity using Athy (1930)
        por = permeability.Permeability.Porosity(rock, depth)
    
        # calculate rock permeability using multipoint method (Hantschel 2009)
        kv, kh = permeability.Permeability.k(rock, por)      
        # calculate connate water saturation using Holmes (2009)
        SW = permeability.Permeability.SwiZ(rock, depth)
            
        # extract values and uncertainties
        Sw = SW[2].n                            # Sw value as float
        Swuc = SW[2].s                          # Sw uncertainty as float
    
        # relative water-liquid permeability with uncertainties
        krel = permeability.Permeability.krp(rock, "WL", ufloat(Sw, Swuc)) 
            
        # extract values and uncertainties
        kr = krel[1].n                          # relative perm water
        kr_uc = krel[1].s                       # relative perm water error
        krow = krel[2].n                        # relative perm oil-water
        krow_uc = krel[2].s                     # relative perm oil-water error
  
        mD = 9.869233e-15                       # conversion from mD to m^2

        ############
        # relative water-liquid permeability with uncertainties
        krell = permeability.Permeability.krp(rock, "VL", ufloat(Sw, Swuc)) 
    
    # extract values and uncertainties
        krg = krell[1].n                         # relative perm water
        krg_uc = krell[1].s                      # relative perm water error
        krog = krell[2].n                        # relative perm oil-water
        krog_uc = krell[2].s                     # relative perm oil-water error   
    
        # Aziz and Settari (1979) approximation - phases do not interact
        Kro = ufloat(krow, krow_uc)*ufloat(krog, krog_uc)
        kro = Kro.n
        kro_uc = Kro.s

        # use oil-water permeability for vertical and horizontal keffs
        Kveff = 10**(kv) * mD * ufloat(kro, kro_uc)
        kveff = Kveff.n
        kveff_uc = Kveff.s    
        Kheff = 10**(kh) * mD * ufloat(kro, kro_uc)
        kheff = Kheff.n
        kheff_uc = Kheff.s
        #############
       
        # calculate temperature and pressure from depth for normal geological
        # conditions
        
        pt = data.Data.PT(depth, tsurf)
        T = pt[0]                               # temperature in kelvin
        P = pt[1]*1e6                           # *1e6 for pressure in Pa

        #print(T,P)  

            # ViscosityPure changed to include mu_l and mu_g on 27/07/21!
        if len(name.MWs) == 1 and name != "H2O":
            visc = viscosity.Viscosity.Pure(name, EOS, Vc, T, P)
            phase = visc[2]
            mul = visc[8]
            mug = visc[9]
            rhol = visc[11]
            rhog = visc[13]
            p = visc[4]
            t = visc[5]
            v_L = visc[6]
            v_G = visc[7]

    
        elif len(name.MWs) == 1 and name == "H2O":
            visc = viscosity.Viscosity.Pure(name, EOS, Vc, T, P)
            phase = visc[2]
            mul = visc[8]
            mug = visc[9]
            rhol = visc[11]
            rhog = visc[13]
            p = visc[4]
            t = visc[5]     
   
        if phase == 'l/g':
            # verticql velocity (liquid)
            mob_vl = 1*(ufloat(kveff, kveff_uc)/mul) # -1 omitted         
            # vertical velocity (gas)
            mob_vg = 1*(ufloat(kveff, kveff_uc)/mug) # -1 omitted      
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vl.n, mob_vl.s, \
                    mob_vg.n, mob_vg.s, rhol, rhog
               
        elif phase == 'l':
            # vertical mobility
            mob_vl = 1*(ufloat(kveff, kveff_uc)/mul) # -1 omitted
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vl.n, mob_vl.s, rhol

        elif phase == 'g':
            # vertical velocity
            mob_vg = 1*(ufloat(kveff, kveff_uc)/mug) # -1 omitted      
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vg.n, mob_vg.s, rhog        

class Run:

    def run(rock, depth, tsurf, eos, plot):

        #depths = []
        mobs = []
        rhows = []
        buoys = []
        vels = []
        dens = []
        viscs = []        

        for z in depth:
    
            mob = Mobility.Mobility("H2", rock, z, tsurf, eos)
            rhow = Mobility.Mobility("H2O", rock, z, tsurf, eos)
            density = mob[-1]
            mobb = mob[8]

            pt = data.Data.PT(z, tsurf)
            T = pt[1]                               # temperature in kelvin
            P = pt[2]*1e6                           # *1e6 for pressure in Pa

            name = data.Fluid.Name(fluid)[0]
            Vc = data.Fluid.Name(fluid)[1]
            mu = viscosity.Viscosity.Pure(name, eos, Vc, T, P)
            if mu[2] == 'g':
                visc = mu[9] 
            elif mu[2] == 'l':
                visc = mu[8]
            elif mu[2] == 'l/g':
                visc = mu[8] # choose liquid viscosity in case of 'l/g' phase
            
            
            
            # correct to ensure liquid density is chosen for  case of l/g phase
            if rhow[5] == 'l/g':
                rhow = rhow[-2]
            elif rhow[5] != 'l/g':
                rhow = rhow[-1]

            buoy = 9.08665*(rhow - mob[-1])
            vel = mob[8]*buoy*3.154e7 # multiply by 3.154e7 s in a year
                
            mobs.append(mobb)
            rhows.append(rhow)
            buoys.append(buoy)
            vels.append(vel)
            dens.append(density)
            viscs.append(visc)

            dict = {'Depth [km]':depth,'vmax [m/year]':vels}
            df = pd.DataFrame(dict)
            print("Mobility algorithm results for H2 and", rock)
            print(df)

        if plot == "on":
        
            # Plotting
            plt.figure(figsize=(8, 6))
            plt.scatter(vels, depth)
            # plot.plot(vels, depth) # plot as line

            plt.ylim(0, max(depth)+ 0.5)  # Depth range from 0 to 4 km
            plt.xlim(0, max(vels) + 2)  # X-axis range from 0 to 2 units greater than the largest velocity
            plt.gca().invert_yaxis()

            
            plt.xlabel('Vertical Velocity [m/year]')
            plt.ylabel('Depth [km]')
            plt.title('Vertical Velocity vs. Depth')
            plt.grid(True)

            plt.show()      

        elif plot == "off":
            pass