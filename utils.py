import numpy as np
import pickle
import os
import warnings
from photochem import Atmosphere, zahnle_earth
from photochem._photochem import PhotoException

class FluxRetrieval():

    mxsteps = 50000
    equilibrium_time = 1e17
    burnin_redox_steps = 5000
    check_redox_steps = 1000
    redox_tolerance = 1e-5
    species = ['CO2','CO','CH4','H2','O2']
    
    def __init__(self, a, b, c, d, output_file):
        self.pc = Atmosphere(a,b,c,d)

        if os.path.exists(output_file):
            warnings.warn(output_file+' already exists. This simulation will append to it.')
        else:
            with open(output_file, 'wb') as f:
                pass
    
        self.output_file = output_file
        
    def equilibrium(self, usol):
        self.pc.initialize_stepper(usol)
        tn = 0.0
        i = 0
        j = 0
        success = True
        try:
            while tn < self.equilibrium_time:
                tn = self.pc.step()
                j += 1
                i += 1
                if i > self.check_redox_steps and j > self.burnin_redox_steps:
                    if self.pc.redox_conservation() < self.redox_tolerance:
                        success = True
                        break
                    else:
                        i = 0
                if j > self.mxsteps:
                    if self.pc.redox_conservation() < self.redox_tolerance:
                        success = True
                        break
                    else:
                        raise PhotoException('too many steps')
        except PhotoException:
            success = False
        self.pc.destroy_stepper()
        
        return success

    def column_average_mix(self):
        # total column
        dz = self.pc.var.z[1] - self.pc.var.z[0]
        col = np.sum(self.pc.wrk.density*dz)

        sol = self.pc.mole_fraction_dict()
        out = {}
        for sp in sol:
            if sp in self.pc.dat.species_names:
                sp_col = np.sum(self.pc.wrk.density*sol[sp]*dz)
                mix = sp_col/col
                out[sp] = mix

        return out
    
    def equilibrium_result(self, success):
        fluxes = self.pc.gas_fluxes()[0]
        sol = self.pc.mole_fraction_dict()
        mix = self.column_average_mix()
        
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

        out['fO2'] = sol['O2'][0]
        out['FO2'] = fluxes['O2']
        out['vO2'] = -fluxes['O2']/(sol['O2'][0]*self.pc.wrk.density[0])

        # column average O3 mixing ratio
        out['cO3'] = mix['O3']

        if not success:
            for key in out:
                out[key] = np.nan

        return out
    
    def save_output(self, success, res):
        out = success, res
        with open(self.output_file,'ab') as f:
            pickle.dump(out, f)
    
    def find_equilibrium(self, log10mix) -> None:
        mix = 10.0**log10mix
        CO2, CO, CH4, H2, O2 = mix
        self.pc.set_lower_bc('CO2',bc_type='mix',mix = CO2)
        self.pc.set_lower_bc('CO',bc_type='mix',mix = CO)
        self.pc.set_lower_bc('CH4',bc_type='mix',mix = CH4)
        self.pc.set_lower_bc('H2',bc_type='mix',mix = H2)
        self.pc.set_lower_bc('O2',bc_type='mix',mix = O2)
        
        # first try the correct answer
        success = self.equilibrium(self.pc.var.usol_init)
        if success:
            self.save_output(success, self.equilibrium_result(success))
            return
        
        # next try unform mixing ratios
        usol = np.ones_like(self.pc.var.usol_init)*1e-40
        for i,sp in enumerate(self.species):
            ind = self.pc.dat.species_names.index(sp)
            usol[ind,:] = mix[i]
        success = self.equilibrium(usol)
        if success:
            self.save_output(success, self.equilibrium_result(success))
            return

        # next try empty atmosphere
        usol = np.ones_like(self.pc.var.usol_init)*1e-40
        success = self.equilibrium(usol)
        self.save_output(success, self.equilibrium_result(success))