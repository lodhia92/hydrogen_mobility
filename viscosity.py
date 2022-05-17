# Calculate viscosity of a chemical mixture using LBC method
# Bhavik Harish Lodhia 08/06/2021

from thermo import ChemicalConstantsPackage
from thermo.eos_mix import PR78MIX, TWUPRMIX, SRKMIX, TWUSRKMIX, APISRKMIX, RKMIX
from thermo.eos_mix import VDWMIX
from chemicals.viscosity import Lorentz_Bray_Clarke, mu_IAPWS
from chemicals.iapws import iapws95_rho
from chemicals import Vm_to_rho
import numpy as np


class Viscosity:
    
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
        mole_fractions = Viscosity.xmole(ws, Mw)
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

    def Pure(name, equation, Vc, t, p):
    
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
    def Mix(name, equation, t, p, Vcs, ws, rhos):
    
        zs = Viscosity.xmole(ws, name.MWs)
        # print(zs)
        
        # Solve EOS to calculate liquid and gas volumes
        if equation == "PR78":
            EOS = PR78MIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                          zs=zs)
            if EOS.phase == str('l'):
                # rhol = Rhol(ws, rhos)    
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))               
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))   
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))              
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)                          
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)
    
        elif equation == "TWUPR":
            zs = Viscosity.xmole(ws, name.MWs)        
            EOS = TWUPRMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                           zs=zs)   
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                     
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))                 
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))               
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)            
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)      
    
        elif equation == "SRK":
            zs = Viscosity.xmole(ws, name.MWs)          
            EOS = SRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                         zs=zs)
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                     
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))                 
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))              
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)               
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)   
                
        elif equation == "TWUSRK":
            zs = Viscosity.xmole(ws, name.MWs)          
            EOS = TWUSRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                            zs=zs)
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                         
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))                  
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))           
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)              
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)   
   
        elif equation == "APISRK":
            zs = Viscosity.xmole(ws, name.MWs)          
            EOS = APISRKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                            zs=zs) 
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                  
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))              
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))               
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)              
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog) 
    
        elif equation == "RK":
            zs = Viscosity.xmole(ws, name.MWs)          
            EOS = RKMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                        zs=zs)
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                   
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))             
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))            
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)  
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)                 
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)  
    
        elif equation == "VDW":
            zs = Viscosity.xmole(ws, name.MWs)          
            EOS = VDWMIX(T=t, P=p, Tcs=name.Tcs, Pcs=name.Pcs, omegas=name.omegas, 
                         zs=zs)
            if EOS.phase == str('l'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))                       
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs) 
                return (name.names, equation, EOS.phase, str('Volume'), p, t,  
                        EOS.V_l, str('NaN'), mu_l, str('NaN'), str('rhol'), rhol,
                        str('rhog'), str('NaN'))
            elif EOS.phase == str('g'):
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))                
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)         
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        str('NaN'), EOS.V_g, str('NaN'), mu_g, str('rhol'), 
                        str('NaN'), str('rhog'), rhog)
            elif EOS.phase == str('l/g'):
                rhol = Vm_to_rho(Vm=EOS.V_l, MW=Viscosity.MwMix(ws, name.MWs))  
                rhog = Vm_to_rho(Vm=EOS.V_g, MW=Viscosity.MwMix(ws, name.MWs))               
                mu_l = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_l, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)
                mu_g = Lorentz_Bray_Clarke(T=t, P=p, Vm=EOS.V_g, zs=zs, MWs=name.MWs, 
                                           Tcs=name.Tcs, Pcs=name.Pcs, Vcs=Vcs)            
                return (name.names, equation, EOS.phase, str('Volume'), p, t, 
                        EOS.V_l, EOS.V_g, mu_l, mu_g, str('rhol'), rhol,
                        str('rhog'), rhog)  
        
    def H2O(equation, t, p):
    
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
      

#name, vc = data.Fluid.Name("CO2")
#equation = "PR78"
#for T in range(300,510,10):
#    for P in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:    
#     # for P in [1, 10, 20, 30, 40]:    
#        Mu = Viscosity.Pure(name, equation, vc, T, P*1e6)
#        print(Mu)
        
#name, vc = data.Fluid.Name("H2O")
#Name = "H2O"
#equation = "PR78"
#for T in range(300,505,5):
#    for P in [1, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:    
#        Mu = Viscosity.H2O(equation, T, P*1e6)
#        print(Name, Mu)

































