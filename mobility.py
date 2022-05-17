# Calculate fluid velocities as a function of rock type, pressure and
# temperature
# By Bhavik Harish Lodhia
# UNSW Sydney 14/07/2021

from viscosity import *             # includes volume functions
from permeability import *          # includes rock and relative permeability


##############################################################################
# Critical volumes , /1000 since we need m^3/mol not m^3/kmol
methane_vc = [0.0986/1000]
CO2_vc = [0.0940/1000]
H2_vc = [0.0642/1000]
H2O_vc = [0.0560/1000]
drygas_vc = [0.0977/1000]
wetgas_vc = [0.1345/1000] 
voil_vc = [0.1956/1000]
lightoil_vc = [0.1895/1000]
blackoil_vc = [0.3121/1000]
mediumoil_vc = [0.3302/1000]
heavyoil_vc = [0.4080/1000]

    
# Define chemicals
methane = ChemicalConstantsPackage(MWs=[16.043], names=['methane'],
          omegas=[0.0115], Pcs=[4.599e6], Tcs=[190.25])
CO2 = ChemicalConstantsPackage(MWs=[44.010], names=['CO2'], omegas=[0.2276], 
      Pcs=[7.382e6], Tcs=[304.19])
H2 = ChemicalConstantsPackage(MWs=[2.016], names=['H2'], omegas=[-0.2150], 
     Pcs=[1.313e6], Tcs=[33.18])
H2O = ChemicalConstantsPackage(MWs=[18.015], names=['H2O'], omegas=[0.3449], 
      Pcs=[22.055e6], Tcs=[647.13])
drygas = ChemicalConstantsPackage(MWs=[17.943], names=['dry_gas'], 
         omegas=[0.0221], Pcs=[4.850e6], Tcs=[197.42])
wetgas = ChemicalConstantsPackage(MWs=[30.186], names=['wet_gas'], 
         omegas=[0.0624], Pcs=[4.801e6], Tcs=[272.40])
H2O = ChemicalConstantsPackage(MWs=[18.015], names=['H2O'], omegas=[0.3449], 
      Pcs=[22.055e6], Tcs=[647.13])


volatileoil = ChemicalConstantsPackage(MWs=[53.135], names=['voloil'],
           omegas=[0.1267],Pcs=[4.373e6], Tcs=[367.79])
lightoil = ChemicalConstantsPackage(MWs=[48.439], names=['lightoil'],
           omegas=[0.1238],Pcs=[4.389e6], Tcs=[357.38])
blackoil = ChemicalConstantsPackage(MWs=[90.072], names=['blackoil'],
           omegas=[0.2423],Pcs=[3.535e6], Tcs=[491.62])
mediumoil = ChemicalConstantsPackage(MWs=[101.141], names=['mediumoil'],
           omegas=[0.2606],Pcs=[3.459e6], Tcs=[511.92])
heavyoil = ChemicalConstantsPackage(MWs=[127.492], names=['heavyoil'],
           omegas=[0.3290],Pcs=[3.041e6], Tcs=[568.42])


mixture = ChemicalConstantsPackage(MWs=[16.043, 30.070, 44.096, 58.123, 
72.150, 84.000, 209.610], names=['methane', 'ethane', 'propane', 'butane',
'pentane', 'C6', 'C7+'], omegas=[0.0115, 0.0995, 0.1523, 0.2002, 0.2515,
0.2510, 0.5295], Pcs=[4.599e6, 4.872e6, 4.248e6, 3.796e6, 3.370e6, 3.217e6,
1.700e6], Tcs=[190.25, 305.32, 369.83, 425.12, 469.70, 510.00, 773.15])
##############################################################################


def fluid(Fluid):

    if Fluid == "CO2":
        name = CO2
        Vc = CO2_vc
        
    elif Fluid == "methane":
        name = methane
        Vc = methane_vc        

    elif Fluid == "H2":
        name = H2
        Vc = H2_vc
        
    elif Fluid == "drygas":
        name = drygas
        Vc = drygas_vc

    elif Fluid == "wetgas":
        name = wetgas
        Vc = wetgas_vc 
        
    elif Fluid == "voil":
        name = volatileoil
        Vc = voil_vc
    
    elif Fluid == "lightoil":
        name = lightoil
        Vc = lightoil_vc
        
    elif Fluid == "blackoil":
        name = blackoil
        Vc = blackoil_vc
    
    elif Fluid == 'mediumoil':
        name = mediumoil
        Vc = mediumoil_vc
        
    elif Fluid == "heavyoil":
        name = heavyoil
        Vc = heavyoil_vc    
        
    elif Fluid == "H2O":
        name = H2O
        Vc = H2O_vc           
        
    # return name, Vc

    elif Fluid == "mixture":
        name = mixture
        Vc = [0.0986/1e3, 0.1455/1e3, 0.2/1e3, 0.255/1e3, 0.313/1e3, 0.3047/1e3, 
              0.7551/1e3]
        ws = [0.5, 0.06, 0.05, 0.03, 0.02, 0.02, 0.32]
        rhos = [300., 356.2, 507., 584., 631., 663.8, 737.3]
        rho = Rhol(ws, rhos)             # density in kg/m^3
        return name, Vc, ws, rhos, rho        
    
    if Fluid == "CO2" or "methane" or "H2" or "drygas" or "wetgas" or "voil" \
                or "lightoil" or "mediumoil" or "blackoil" or "heavyoil":                        
        return name, Vc
        
        
def buoyancy(Fluid, EOS, depth, tsurf):
    
    f = fluid(Fluid)
    name = f[0]
    Vc = f[1]
    
    if len(name.MWs) > 1:
        # f = fluid(Fluid)
        # name = f[0]
        # Vc = f[1]
        ws = f[2]
        rhos = f[3]
        rho = f[4]

    # calculate temperature and pressure from depth for normal geological
    # conditions

    pt = PT(3, depth, tsurf)
    T = pt[1]                               # temperature in kelvin
    P = pt[2]*1e6                           # *1e6 for pressure in Pa
    # return depth, T, P
        
    # ViscosityMix needs updating to correct mu_l and mu_g!
    if len(name.MWs) > 1:
        visc = ViscosityMix(name, EOS, T, P, Vc, ws, rhos)
        phase = visc[2]
        p = visc[4]
        t = visc[5]        
        mul = visc[8]
        mug = visc[9]
        rhol = visc[11]
        rhog = visc[13] 
        # print(phase, rhol, rhog)
        
    # ViscosityPure changed to include mu_l and mu_g on 27/07/21!
    elif len(name.MWs) == 1:
        visc = ViscosityPure(name, EOS, Vc, T, P)
        phase = visc[2]
        p = visc[4]
        t = visc[5]       
        mul = visc[8]
        mug = visc[9]
        rhol = visc[11]
        rhog = visc[13]
        # print(phase, rhol, rhog)
    # return depth, T, P, mu, rhol, rhog
    
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
def mobility(Fluid, rock, depth, tsurf, EOS):
    
    f = fluid(Fluid)
    name = f[0]
    Vc = f[1]
    
    if len(name.MWs) > 1:
        # f = fluid(Fluid)
        # name = f[0]
        # Vc = f[1]
        ws = f[2]
        rhos = f[3]
        rho = f[4]
        
    # calculate porosity using Athy (1930)
    por = porosity(rock, depth)             # porosity from Athy (1930)
    # return por
    
    # calculate rock permeability using multipoint method (Hantschel 2009)
    krock = k("multipoint", rock, por)      # rock permeability
    # return krock
    kv = krock[0]                           # vertical rock permeability
    kh = krock[1]                           # horizontal rock permeability
    # return kv, kh
    # calculate connate water saturation using Holmes (2009)
    SW = SwiZ(rock, depth)                     # connate water saturation
    
    # extract values and uncertainties
    Sw = SW[2].n                            # Sw value as float
    Swuc = SW[2].s                          # Sw uncertainty as float
    # return Sw, Swuc
    
    # relative water-liquid permeability with uncertainties
    krel = krp(rock, "Quadratic", "WL", ufloat(Sw, Swuc)) 
            
    # extract values and uncertainties
    kr = krel[1].n                          # relative perm water
    kr_uc = krel[1].s                       # relative perm water error
    krow = krel[2].n                        # relative perm oil-water
    krow_uc = krel[2].s                     # relative perm oil-water error
    # return kr, kr_uc, krow, krow_uc
  
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
    krell = krp(rock, "Quadratic", "VL", ufloat(Sw, Swuc)) 
    
    # extract values and uncertainties
    krg = krell[1].n                         # relative perm water
    krg_uc = krell[1].s                      # relative perm water error
    krog = krell[2].n                        # relative perm oil-water
    krog_uc = krell[2].s                     # relative perm oil-water error
    # return krg, krg_uc, krog, krog_uc    
    
    # Aziz and Settari (1979) approximation - phases do not interact
    Kro = ufloat(krow, krow_uc)*ufloat(krog, krog_uc)
    kro = Kro.n
    kro_uc = Kro.s
    # return kro, kro_uc

    # use oil-water permeability for vertical and horizontal keffs
    Kveff = 10**(kv) * mD * ufloat(kro, kro_uc)
    kveff = Kveff.n
    kveff_uc = Kveff.s    
    Kheff = 10**(kh) * mD * ufloat(kro, kro_uc)
    kheff = Kheff.n
    kheff_uc = Kheff.s
    # return kveff, kveff_uc, kheff, kheff_uc
    #############
       
    # calculate temperature and pressure from depth for normal geological
    # conditions

    pt = PT(3, depth, tsurf)
    T = pt[1]                               # temperature in kelvin
    P = pt[2]*1e6                           # *1e6 for pressure in Pa
    # return depth, T, P

    # define variables for inputs    
       
    # calculate temperature and pressure from depth for normal geological
    # conditions

    pt = PT(3, depth, tsurf)
    T = pt[1]                               # temperature in kelvin
    P = pt[2]*1e6                           # *1e6 for pressure in Pa
    # return depth, T, P
        
    if len(name.MWs) > 1:
        visc = ViscosityMix(name, EOS, T, P, Vc, ws, rhos)
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
        
        # if phase == 'l':
        #     rhol = visc[10]
        # elif phase == 'g':
        #     rhog = visc[10]
        # elif phase == 'l/g':
        #     rhol = visc[8]
        #     rhog = visc[10]        
        
    # ViscosityPure changed to include mu_l and mu_g on 27/07/21!
    elif len(name.MWs) == 1 and fluid != H2O:
        visc = ViscosityPure(name, EOS, Vc, T, P)
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
    
    elif len(name.MWs) == 1 and fluid == H2O:
        visc = ViscosityPure(name, EOS, Vc, T, P)
        phase = visc[2]
        mul = visc[8]
        mug = visc[9]
        rhol = visc[11]
        rhog = visc[13]
        p = visc[4]
        t = visc[5]     
    
    g = 9.08665                             # gravity in m/s^2 
    rhow = iapws95_rho(T=t, P=p)   
   
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

#######




depths = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6, \
          1.7, 1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2, \
          3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0]

# equations = ["PR78", "SRK", "RK"]
# tsurfs = [0., 20.]
# for tsurf in tsurfs:   
#     for EOS in equations:
#         for z in depths:
#             out = buoyancy("CO2", EOS, z, tsurf)
#             print(out)
            # out = buoyancy("methane", EOS, z, tsurf)
            # print(out)   
            # out = buoyancy("H2", EOS, z, tsurf)
            # print(out)    
            # out = buoyancy("mixture", EOS, z, tsurf)
            # print(out)
            # out = buoyancy("drygas", EOS, z, tsurf)
            # print(out)   
            # out = buoyancy("wetgas", EOS, z, tsurf)
            # print(out)    
            # out = buoyancy("voil", EOS, z, tsurf)
            # print(out)   
            # out = buoyancy("lightoil", EOS, z, tsurf)
            # print(out)    
            # out = buoyancy("blackoil", EOS, z, tsurf)
            # print(out)   

rocks = ["Sandstone"] 
# rocks = ["Micrite"]           
equations = ["PR78", "SRK", "RK"]
HCs = ["PR78", "SRK", "TWUPR", "TWUSRK"]
tsurfs = [0., 20.]

for tsurf in tsurfs:   
    for EOS in equations:
        for z in depths:
            for rock in rocks:
                out = mobility("CO2", rock, z, tsurf, EOS)
                print(out)
                # out = mobility("methane", rock, z, tsurf, EOS)
                # print(out)
                # out = mobility("H2", rock, z, tsurf, EOS)
                # print(out)      
                # out = mobility("drygas", rock, z, tsurf, EOS)
                # print(out)  
                # out = mobility("wetgas", rock, z, tsurf, EOS)
                # print(out)                
                

# for tsurf in tsurfs:   
#     for EOS in HCs:
#         for z in depths:
#             for rock in rocks:
                # out = mobility("mixture", rock, z, tsurf, EOS)
                # print(out)                    
                # out = mobility("voil", rock, z, tsurf, EOS)
                # print(out)                
                # out = mobility("lightoil", rock, z, tsurf, EOS)
                # print(out) 
                # out = mobility("blackoil", rock, z, tsurf, EOS)
                # print(out)    
            
            
