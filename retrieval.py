import numpy as np
import pickle
from multiprocessing import Pool
import dynesty
import rfast

# initialize radiative transfer and retrieval
r = rfast.Rfast('input/inputs.scr')
r.initialize_retrieval("input/rpars.txt")

# Make fake data of Modern Earth twin.
F1, F2 = r.genspec_scr()
dat, err = r.noise(F2)

# functions that need to be fed to nested sampler.
# we pass data as globals.
def lnlike(x_t):
    return rfast.lnlike_nest(x_t, r, dat, err)

def prior_transform(u):
    return rfast.prior_transform(u, r)
        
if __name__ == "__main__":
    # number of processes
    NPROCESS = 4
    
    # retrieval with all gases
    with Pool(NPROCESS) as pool:
        sampler = dynesty.NestedSampler(lnlike, prior_transform, r.retrieval.nret, \
                                        nlive = 1000, pool=pool, queue_size=NPROCESS)
        sampler.run_nested(dlogz=0.001)

    with open('normal_retrieval.pkl', 'wb') as f:
        pickle.dump(sampler.results, f)
