
```sh
# create conda environment
conda create -n rfast -c conda-forge python numpy numba scipy ruamel.yaml astropy emcee dynesty multiprocess scikit-build cython pyyaml threadpoolctl p-tqdm

# activate
conda activate rfast

# rfast
git clone https://github.com/Nicholaswogan/rfast.git
cd rfast
git checkout b9e5bc8cb27770448fd998b0ce66ebdf13ec4ea3
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf rfast

# photochem
git clone --recursive --branch=dev https://github.com/Nicholaswogan/photochem
cd photochem
git checkout 3cb0cdf636dc1dddc9934d322c7165edc115be71
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf photochem

# clima
git clone --recursive --branch=dev https://github.com/Nicholaswogan/clima
cd clima
git checkout c6b05bdf2d5eccdc9db124e0d3fe6d04d8e44c55
python -m pip install --no-deps --no-build-isolation . -v
cd ..
rm -rf clima
```