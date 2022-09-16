import numpy as np
import pickle
from dynesty import utils as dyfunc

import rfast
from photochem import zahnle_earth
from utils import FluxRetrieval
from p_tqdm import p_umap
from threadpoolctl import threadpool_limits
threadpool_limits(limits=1)

np.random.seed(0)

def get_samples(results, nr):
    samples = results['samples']
    weights = np.exp(results['logwt'] - results['logz'][-1])
    new_samples = dyfunc.resample_equal(samples, weights)
    inds = np.random.randint(low=0, high=new_samples.shape[0]-1, size=nr)
    samples = np.empty((nr,new_samples.shape[1]))
    for i in range(new_samples.shape[1]):
        samples[:,i] = new_samples[inds,i]
    return new_samples, samples

def samples_for_photo(filename, nr, flx, r):
    with open(filename,'rb') as f:
        res = pickle.load(f)
    all_samp, samp = get_samples(res, nr)
    samp_n = np.empty((nr,len(flx.species)))
    for i,sp in enumerate(flx.species):
        ind = list(r.retrieval.param_names).index('f'+sp.lower())
        samp_n[:,i] = samp[:,ind]
    return [samp_n[i,:] for i in range(samp_n.shape[0])]

# stuff to be passed by scope
r = rfast.Rfast('input/inputs_Proterozoic.scr')
r.initialize_retrieval("input/rpars.txt")
flx = FluxRetrieval(zahnle_earth,\
            "./input/settings_Proterozoic.yaml",\
            "./input/Sun_2.4Ga.txt",\
            "./input/atmosphere_Proterozoic.txt",\
            "results/Proterozoic_flux_retrieval_1.pkl")
flx.mxsteps = 30000
flx.pc.var.verbose = 0

def wrapper(log10mix):
        return flx.find_equilibrium(log10mix)

if __name__ == "__main__":
    filename = 'results/Proterozoic_retrieval_1.pkl'
    nr = 100000
    NUM_PROCESS = 48

    samp = samples_for_photo(filename, nr, flx, r)
    p_umap(wrapper, samp, num_cpus=NUM_PROCESS)

