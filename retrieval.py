import numpy as np
import pickle
from multiprocessing import Pool
import dynesty
import rfast

# initialize radiative transfer and retrieval
r = rfast.Rfast('input/inputs_Proterozoic_1.scr')
r.initialize_retrieval("input/rpars.txt")

F1, F2 = r.genspec_scr()
dat, err = r.noise(F2)

def lnlike(x_t):
    return rfast.lnlike_nest(x_t, r, dat, err)

def prior_transform(u):
    return rfast.prior_transform(u, r)
        
if __name__ == "__main__":
    # number of processes
    NPROCESS = 48
    
    # retrieval with all gases
    with Pool(NPROCESS) as pool:
        sampler = dynesty.NestedSampler(lnlike, prior_transform, r.retrieval.nret, \
                                        nlive = 1000, pool=pool, queue_size=NPROCESS)
        sampler.run_nested(dlogz=0.001)

    with open('results/Proterozoic_retrieval_1.pkl', 'wb') as f:
        pickle.dump(sampler.results, f)
