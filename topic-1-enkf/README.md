# [Ensemble Kalman Filter](http://lumetta.web.engr.illinois.edu/408-S20/projects/)

Libaries and dependencies are installed in the docker image already. They are installed under `$ECE408_BASE` (i.e., `/ece408`) inside the container.

See [instructions](http://lumetta.web.engr.illinois.edu/408-S20/projects/Kalman-filter.pdf) and the [Dockerfile](./Dockerfile) for detailed make/installation process if you plan to develop locally.

## Implementation
As mentioned in the [instruction](http://lumetta.web.engr.illinois.edu/408-S20/projects/Kalman-filter.pdf), you'll accelerate the core computation that is equvilent to the following in python/numpy.
```python
P_HT = np.multiply(C, e @ e.T) @ H.T / (L - 1)
```

#### About `src`
Several python scripts are included in [`lenkf/python`](lenkf/python). They are also copied to [`src/`](src) (with minor modification in [`compute_P_HT.py`](src/compute_P_HT.py) to set the correct path). The script generates the input matrices (to `/tmp`) and compares result of the C version against the numpy version for `compute_P_HT`. `compute_P_HT` is built from [`lenkf`](lenkf) which uses [`libmdb_matrix`](libmdb_matrix). Check [`lenkf.c`](lenkf/lenkf.c) (and other source code) to see detailed implementation. Feel free to change or reuse code snippets in the libraries for boilerplate in your CUDA-accelerated version.

Under [`src/`](src) directory, you can find a RAI configuration file `rai_build.yml`. Currently it runs `compute_P_HT.py` and copy all the generated matrices to the current working directory. Make changes to the commands as you start working and building the CUDA-acceleratred version.
