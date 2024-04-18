# Calculate viscosity of a chemical mixture using LBC method
# Bhavik Harish Lodhia

from thermo import ChemicalConstantsPackage
from thermo.eos_mix import PR78MIX, TWUPRMIX, SRKMIX, TWUSRKMIX, APISRKMIX, RKMIX
from thermo.eos_mix import VDWMIX
from chemicals.viscosity import Lorentz_Bray_Clarke, mu_IAPWS
from chemicals.iapws import iapws95_rho
from chemicals import Vm_to_rho
import numpy as np


class Viscosity:

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
































