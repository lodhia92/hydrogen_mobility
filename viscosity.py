# Calculate viscosity of a chemical mixture using LBC method
# Bhavik Harish Lodhia 08/06/2021

from thermo import ChemicalConstantsPackage
from thermo.eos_mix import PR78MIX, TWUPRMIX, SRKMIX, TWUSRKMIX, APISRKMIX, RKMIX
from thermo.eos_mix import VDWMIX
from chemicals.viscosity import Lorentz_Bray_Clarke, mu_IAPWS
from chemicals.iapws import iapws95_rho, iapws97_rho
from chemicals import Vm_to_rho
import numpy as np

# Critical volumes 
methane_vc = [0.0986/1000]
CO2_vc = [0.0940/1000]
H2_vc = [0.0642/1000]
H2O_vc = [0.0560/1000]
drygas_vc = [0.0977/1000]
wetgas_vc = [0.1345/1000] 
voil_vc = [0.1956/1000]
lightoil_vc = [0.1895/1000]
backoil_vc = [0.3121/1000]
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

# # Calculate mole fraction for Hantschel2009 from weight fractions
# ws = [0.5, 0.06, 0.05, 0.03, 0.02, 0.02, 0.32]       # mass fractions of components
# ni = [m*n for m, n in zip(H2009.MWs, ws)]            # calulate moles of components
# ntot = [sum(ni)]*len(ws)                             # calculate total moles
# H2009_zs = [a/b for a,b in zip(ni, ntot)]            # calculate molar fraction

#     [methane, ethane, propane, butane, pentane, C6, C7+]
mixture_vcs = [0.0986/1e3, 0.1455/1e3, 0.2/1e3, 0.255/1e3, 0.313/1e3, 0.3047/1e3, 
            0.7551/1e3]


def xmole(ws, Mw):
    
    if len(ws) != len(Mw):
        print('Weight fractions and molecular weights must have same number \
               of components!')
        pass
    # elif sum(ws) != 1.:
    #     print('Sum of mass fractions not equal to 1 but is ', sum(ws))
    #     pass

    wx = np.array(ws)*100   # multiply by 100 for use with molar equation
    Mws = np.array(Mw)      # array of molecular masses

    n = wx/Mws              # moles of each component
    ntot = sum(n)           # total number of moles
    Zs = n/ntot             # mole fractions (array)
    # print(Zs)
    
    zs = Zs.tolist()        # return mole fraction as list for chemicals
    return zs
    # print(zs, sum(zs))
    
def MwMix(ws, Mw):
    mole_fractions = xmole(ws, Mw)
    recip = [i / j for i, j in zip(mole_fractions, Mw)]
    Mw_mix = 1./sum(recip)
    return Mw_mix
 
def Rhol(mass_fraction, rho):

    # Create 1D arrays from input lists
    ws = np.array(mass_fraction)
    rhos = np.array(rho)
    
    # print(len(ws))
    
    # Check for input errors
    if len(ws) != len(rhos):      
        print(str('Incorrect number of component properties!'))
        pass
    
    recip = ws/rho  # 1/density = mass fraction / density for each component
    density = 1/sum(recip)
    return density


def EOS_Pure(name, equation, t, p):
    
    # Solve EOS to calculate liquid and gas volumes
    if equation == "PR78":
        eos = PR78MIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        
    elif equation == "TWUPR":
        eos = TWUPRMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])                       
    elif equation == "SRK":
        eos = SRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])       
    elif equation == "TWUSRK":
        eos = TWUSRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
    elif equation == "APISRK":
        eos = APISRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0]) 
    elif equation == "RK":
        eos = RKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
    elif equation == "VDW":
        eos = VDWMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])               

    if eos.phase == str('l'):
        return name.names, equation, eos.phase, str('Volume'), p, t, eos.V_l, str('NaN')
    elif eos.phase == str('l/g'):
        return name.names, equation, eos.phase, str('Volume'), p, t, eos.V_l, eos.V_g
        # print(str(name.names), EOS, str('Volume'), eos.phase, p, t, eos.V_l, eos.V_g)
    elif eos.phase == str('g'):
        return name.names, equation, eos.phase, str('Volume'), p, t, str('NaN'), eos.V_g
        # print(str(name.names), EOS, str('Volume'), eos.phase, p, t, str('NaN'), eos.V_g)
    else:
        pass
    
def EOS_Mix(name, equation, t, p, zs):

    # Solve EOS to calculate liquid and gas volumes
    if equation == "PR78":
        eos = PR78MIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
    elif equation == "TWUPR":
        eos = TWUPRMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)                       
    elif equation == "SRK":
        eos = SRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)       
    elif equation == "TWUSRK":
        eos = TWUSRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
    elif equation == "APISRK":
        eos = APISRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs) 
    elif equation == "RK":
        eos = RKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
    elif equation == "VDW":
        eos = VDWMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)               

    if eos.phase == str('l'):
        return name.names, equation, eos.phase, str('Volume'), p, t, eos.V_l, str('NaN')
    elif eos.phase == str('l/g'):
        return name.names, equation, eos.phase, str('Volume'), p, t, eos.V_l, eos.V_g
        # print(str(name.names), EOS, str('Volume'), eos.phase, p, t, eos.V_l, eos.V_g)
    elif eos.phase == str('g'):
        return name.names, equation, eos.phase, str('Volume'), p, t, str('NaN'), eos.V_g
        # print(str(name.names), EOS, str('Volume'), eos.phase, p, t, str('NaN'), eos.V_g)
    else:
        pass    



def ViscosityPure(name, equation, Vc, t, p):
    
    # Solve EOS to calculate liquid and gas volumes
    if equation == "PR78":
        EOS = PR78MIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)     
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)
    
    elif equation == "TWUPR":
        EOS = TWUPRMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])   
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)           
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)               
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)           
    
    elif equation == "SRK":
        EOS = SRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)              
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)    
    
    elif equation == "TWUSRK":
        EOS = TWUSRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)              
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)   
   
    elif equation == "APISRK":
        EOS = APISRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0]) 
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)               
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)   
    
    elif equation == "RK":
        EOS = RKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)               
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)   
    
    elif equation == "VDW":
        EOS = VDWMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=[1.0])
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                  Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=sum(name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=sum(name.MWs))
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=[1.0], MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vc)               
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)

# Needs updating to correct density bug and mu_l and mu_g for l/g phase    
def ViscosityMix(name, equation, t, p, Vcs, ws, rhos):
    
    zs = xmole(ws, name.MWs)
    # print(zs)
    
    # Solve EOS to calculate liquid and gas volumes
    if equation == "PR78":
        EOS = PR78MIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
        if EOS.phase == str('l'):
            # rhol = Rhol(ws, rhos)    
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))               
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))   
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))              
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)                          
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)
    
    elif equation == "TWUPR":
        zs = xmole(ws, name.MWs)        
        EOS = TWUPRMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)   
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                     
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))                 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))               
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)            
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)      
    
    elif equation == "SRK":
        zs = xmole(ws, name.MWs)          
        EOS = SRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                     
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))                 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))              
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)               
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)   
    
    elif equation == "TWUSRK":
        zs = xmole(ws, name.MWs)          
        EOS = TWUSRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                         
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))                  
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))           
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)              
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)   
   
    elif equation == "APISRK":
        zs = xmole(ws, name.MWs)          
        EOS = APISRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs) 
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                  
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))              
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))               
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)              
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog) 
    
    elif equation == "RK":
        zs = xmole(ws, name.MWs)          
        EOS = RKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                   
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))             
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))            
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)  
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)                 
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)  
    
    elif equation == "VDW":
        zs = xmole(ws, name.MWs)          
        EOS = VDWMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
              zs=zs)
        if EOS.phase == str('l'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))                       
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
            return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                    EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                    str('rhog'), str('NaN'))
        elif EOS.phase == str('g'):
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))                
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                 Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                    str('NaN'), str('rhog'), rhog)
        elif EOS.phase == str('l/g'):
            rhol = Vm_to_rho(Vm=EOS.V_l, MW=MwMix(ws, name.MWs))  
            rhog = Vm_to_rho(Vm=EOS.V_g, MW=MwMix(ws, name.MWs))               
            mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
            mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)            
            return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                    EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                    str('rhog'), rhog)  
        
def ViscosityH2O(equation, t, p):
    
    H2O = ChemicalConstantsPackage(MWs=[18.015], names=['H2O'], omegas=[0.3449], 
          Pcs=[22.055e6], Tcs=[647.13]) 
    # rho = iapws95_rho(T=t, P=p)
    
    # Solve EOS to calculate liquid and gas volumes
    if equation == "PR78":
        EOS = PR78MIX(T=t, P=p, Tcs=H2O.Tcs, Pcs=H2O.Pcs, omegas=H2O.omegas, 
              zs=[1.0])
        rhow = iapws95_rho(T=t, P=p)           
        # rhol = iapws95_rho(T=t, P=p)        
        if EOS.phase == str('l'):               
            mu = mu_IAPWS(T=t, rho=rhow)
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, str('NaN'), 
                    mu, rhow)
        elif EOS.phase == str('g'):         
            mu = mu_IAPWS(T=t, rho=rhow)        
            return (equation, EOS.phase, str('Volume'), p, t, str('NaN'), EOS.V_g, 
                    mu, rhow)
        elif EOS.phase == str('l/g'):            
            mu = mu_IAPWS(T=t, rho=rhow)              
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, EOS.V_g, 
                    mu, rhow)     
    
    elif equation == "SRK":
        EOS = SRKMIX(T=t, P=p, Tcs=H2O.Tcs, Pcs=H2O.Pcs, omegas=H2O.omegas, 
              zs=[1.0])
        rhow = iapws95_rho(T=t, P=p)           
        # rhol = iapws95_rho(T=t, P=p)        
        if EOS.phase == str('l'):               
            mu = mu_IAPWS(T=t, rho=rhow)
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, str('NaN'), 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('g'):            
            mu = mu_IAPWS(T=t, rho=rhow)        
            return (equation, EOS.phase, str('Volume'), p, t, str('NaN'), EOS.V_g, 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('l/g'):            
            mu = mu_IAPWS(T=t, rho=rhow)              
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, EOS.V_g, 
                    mu, rhow)             
    
    elif equation == "RK":
        EOS = RKMIX(T=t, P=p, Tcs=H2O.Tcs, Pcs=H2O.Pcs, omegas=H2O.omegas, 
              zs=[1.0])    
        rhow = iapws95_rho(T=t, P=p)           
        # rhol = iapws95_rho(T=t, P=p)        
        if EOS.phase == str('l'):               
            mu = mu_IAPWS(T=t, rho=rhow)
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, str('NaN'), 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('g'):          
            mu = mu_IAPWS(T=t, rho=rhow)        
            return (equation, EOS.phase, str('Volume'), p, t, str('NaN'), EOS.V_g, 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('l/g'):            
            mu = mu_IAPWS(T=t, rho=rhow)              
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, EOS.V_g, 
                    mu, rhow)   
    
    elif equation == "VDW":
        EOS = VDWMIX(T=t, P=p, Tcs=H2O.Tcs, Pcs=H2O.Pcs, omegas=H2O.omegas, 
              zs=[1.0])
        rhow = iapws95_rho(T=t, P=p)           
        # rhol = iapws95_rho(T=t, P=p)        
        if EOS.phase == str('l'):               
            mu = mu_IAPWS(T=t, rho=rhow)
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, str('NaN'), 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('g'):             
            mu = mu_IAPWS(T=t, rho=rhow)        
            return (equation, EOS.phase, str('Volume'), p, t, str('NaN'), EOS.V_g, 
                    str('Viscosity'), mu, rhow)
        elif EOS.phase == str('l/g'):            
            mu = mu_IAPWS(T=t, rho=rhow)              
            return (equation, EOS.phase, str('Volume'), p, t, EOS.V_l, EOS.V_g, 
                    mu, rhow)      
    
    elif equation == "TWUPR" or equation == "TWUSRK" or equation == "APISRK":
        print(str('Must use PR78, SRK, RK or VDW EOS for water'))      


# name = mixture
# equation = "PR78"
# # vcs = H2009_vcs
# ws = [0.5, 0.06, 0.05, 0.03, 0.02, 0.02, 0.32]
# rhos = [300.0, 356.2, 507.0, 584.0, 631.1, 663.8, 738.4]
# vcs = [0.0986/1e3, 0.1455/1e3, 0.2/1e3, 0.255/1e3, 0.313/1e3, 0.3047/1e3, 0.7551/1e3]
# for T in range(300,505,5):
#     # for P in [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:    
#     for P in [1, 10, 20, 30, 40]: 
#         Mu = ViscosityMix(name, equation, T, P*1e6, vcs, ws, rhos)
#         print(Mu)
      

# name = CO2
# equation = "PR78"
# vcs = CO2_vc
# for T in range(300,510,10):
#     for P in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:    
#     # for P in [1, 10, 20, 30, 40]:    
#         Mu = ViscosityPure(name, equation, vcs, T, P*1e6)
#         print(Mu)
        
# name = H2O
# Name = "H2O"
# equation = "PR78"
# vcs = methane_vc
# for T in range(300,505,5):
#     for P in [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:    
#         Mu = ViscosityH2O(equation, T, P*1e6)
#         print(Name, Mu)

































