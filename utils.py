import numpy as np
from photochem import Atmosphere, zahnle_earth
from photochem._photochem import PhotoException

class FluxRetrieval():
    
    mxsteps = 50000
    equilibrium_time = 1e17
    species = ['CO2','CO','CH4','H2']
    
    def __init__(self,a,b,c,d):
        self.pc = Atmosphere(a,b,c,d)
        
    def equilibrium(self, usol):
        self.pc.initialize_stepper(usol)
        tn = 0.0
        j = 0
        success = True
        try:
            while tn < self.equilibrium_time:
                tn = self.pc.step()
                j += 1
                if j > self.mxsteps:
                    raise PhotoException('too many steps')
        except PhotoException:
            success = False
        self.pc.destroy_stepper()
        
        return success
    
    def equilibrium_result(self, success):
        fluxes = self.pc.gas_fluxes()[0]
        sol = self.pc.mole_fraction_dict()
        
        out = {}
        out['fCO2'] = sol['CO2'][0]
        out['FCO2'] = fluxes['CO2']
        out['vCO2'] = -fluxes['CO2']/(sol['CO2'][0]*self.pc.wrk.density[0])
        
        out['fCO'] = sol['CO'][0]
        out['FCO'] = fluxes['CO']
        out['vCO'] = -fluxes['CO']/(sol['CO'][0]*self.pc.wrk.density[0])
        
        out['fCH4'] = sol['CH4'][0]
        out['FCH4'] = fluxes['CH4']
        out['vCH4'] = -fluxes['CH4']/(sol['CH4'][0]*self.pc.wrk.density[0])
        
        out['fH2'] = sol['H2'][0]
        out['FH2'] = fluxes['H2']
        out['vH2'] = -fluxes['H2']/(sol['H2'][0]*self.pc.wrk.density[0])

        if not success:
            for key in out:
                out[key] = np.nan

        return out
        
    def find_equilibrium(self, log10mix):
        mix = 10.0**log10mix
        CO2, CO, CH4, H2 = mix
        self.pc.set_lower_bc('CO2',bc_type='mix',mix = CO2)
        self.pc.set_lower_bc('CO',bc_type='mix',mix = CO)
        self.pc.set_lower_bc('CH4',bc_type='mix',mix = CH4)
        self.pc.set_lower_bc('H2',bc_type='mix',mix = H2)
        
        # first try the correct answer
        success = self.equilibrium(self.pc.var.usol_init)
        if success:
            return success, self.equilibrium_result(success)
        
        # next try unform mixing ratios
        usol = np.ones_like(self.pc.var.usol_init)*1e-40
        for i,sp in enumerate(self.species):
            ind = self.pc.dat.species_names.index(sp)
            usol[ind,:] = mix[i]
        success = self.equilibrium(usol)
        if success:
            return success, self.equilibrium_result(success)

        # next try empty atmosphere
        usol = np.ones_like(self.pc.var.usol_init)*1e-40
        success = self.equilibrium(usol)
        if success:
            return success, self.equilibrium_result()
        else:
            return success, self.equilibrium_result(success)