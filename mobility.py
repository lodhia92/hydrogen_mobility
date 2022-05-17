# Calculate fluid velocities as a function of rock type, pressure and
# temperature
# By Bhavik Harish Lodhia
# UNSW Sydney

import data
import permeability
import viscosity
from uncertainties import ufloat
from chemicals.iapws import iapws95_rho

class Mobility:

    def Buoyancy(Fluid, EOS, depth, tsurf):
    
        name, Vc = data.Data.Name(Fluid)
    
        if len(name.MWs) == 1:
            visc = viscosity.Viscosity.Pure(name, EOS, Vc, T, P)
            phase = visc[2]
            p = visc[4]
            t = visc[5]       
            mul = visc[8]
            mug = visc[9]
            rhol = visc[11]
            rhog = visc[13]
        
        g = 9.08665                             # gravity in m/s^2 
   
        # calculate water density at PT
        rhow = iapws95_rho(T=t, P=p)   
   
        if phase == 'l/g':
            # liquid buoyancy
            buoy_l = g*(rhow - rhol)
            # gas buoyancy
            buoy_g = g*(rhow - rhog)      
            return Fluid, EOS, tsurf, phase, depth, buoy_l, buoy_g
               
        elif phase == 'l':
            # liquid buoyancy
            buoy_l = g*(rhow - rhol) 
            return Fluid, EOS, tsurf, phase, depth, buoy_l

        elif phase == 'g':
            # gas buoyancy
            buoy_g = g*(rhow - rhog)    
            return Fluid, EOS, tsurf, phase, depth, buoy_g
        
    

    # Calculate fluid mobility as a function of depth
    def Mobility(Fluid, rock, depth, tsurf, EOS):
    
        name, Vc = data.Fluid.Name(Fluid)
            
        # calculate porosity using Athy (1930)
        por = permeability.Permeability.Porosity(rock, depth)
    
        # calculate rock permeability using multipoint method (Hantschel 2009)
        kv, kh = permeability.Permeability.k("multipoint", rock, por)      
        # calculate connate water saturation using Holmes (2009)
        SW = permeability.Permeability.SwiZ(rock, depth)
            
        # extract values and uncertainties
        Sw = SW[2].n                            # Sw value as float
        Swuc = SW[2].s                          # Sw uncertainty as float
    
        # relative water-liquid permeability with uncertainties
        krel = permeability.Permeability.krp(rock, "Quadratic", "WL", ufloat(Sw, Swuc)) 
            
        # extract values and uncertainties
        kr = krel[1].n                          # relative perm water
        kr_uc = krel[1].s                       # relative perm water error
        krow = krel[2].n                        # relative perm oil-water
        krow_uc = krel[2].s                     # relative perm oil-water error
  
        mD = 9.869233e-15                       # conversion from mD to m^2
    
        # Kveff = 10**(kv) * mD * ufloat(kr, kr_uc)
        # kveff = Kveff.n
        # kveff_uc = Kveff.s
        # Kheff = 10**(kv) * mD * ufloat(kr, kr_uc)
        # kheff = Kheff.n
        # kheff_uc = Kheff.s    
        # return kveff, kveff_uc, kheff, kheff_uc
    
        ############
        # relative water-liquid permeability with uncertainties
        krell = permeability.Permeability.krp(rock, "Quadratic", "VL", ufloat(Sw, Swuc)) 
    
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
        
        pt = data.Data.PT(3, depth, tsurf)
        T = pt[1]                               # temperature in kelvin
        P = pt[2]*1e6                           # *1e6 for pressure in Pa
            
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
            # print(phase, rhol, rhog)
            # return depth, T, P, mu, rhol, rhog
    
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
            # horizontal velocity (liquid)
            mob_hl = 1*(ufloat(kheff, kheff_uc)/mul) # -1 omitted          
            # vertical velocity (gas)
            mob_vg = 1*(ufloat(kveff, kveff_uc)/mug) # -1 omitted      
            # horizontal velocity (gas)
            mob_hg = 1*(ufloat(kheff, kheff_uc)/mug) # -1 omitted
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vl.n, mob_vl.s, \
                mob_hl.n, mob_hl.s, mob_vg.n, mob_vg.s, mob_hg.n, mob_hg.s, \
                mul, mug, v_L, v_G, rhol, rhog
               
        elif phase == 'l':
            # vertical mobility
            mob_vl = 1*(ufloat(kveff, kveff_uc)/mul) # -1 omitted
            #horizontal velocity
            mob_hl = 1*(ufloat(kheff, kheff_uc)/mul) # -1 omitted      
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vl.n, mob_vl.s, \
                mob_hl.n, mob_hl.s, mul, v_L, rhol

        elif phase == 'g':
            # vertical velocity
            mob_vg = 1*(ufloat(kveff, kveff_uc)/mug) # -1 omitted    
            # horizontal velocity
            mob_hg = 1*(ufloat(kheff, kheff_uc)/mug) # -1 omitted    
            return Fluid, rock, EOS, depth, tsurf, phase, p, t, mob_vg.n, mob_vg.s, \
                mob_hg.n, mob_hg.s, mug, v_G, rhog      
  
#            
            
