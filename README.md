## Genome-scale metabolic reconstruction (GEM) of *Chlorella ohadii*

Here, we present a fully automated platform for the *de novo* generation of genome-scale metabolic models (GEMs). We deployed this platform to reconstruct the GEM and the enzyme-constrained GEM for *C. ohadii*.

### Contents
Here, the project's directory structure is illustrated: 
```
c_ohadii_GEM
├── code
│   ├── functions      # contains the implementation of the main functions
|   ├── scripts        # includes the python scripts
│   └── examples.m     # provides examples of how to use the functions
├── data               # the data, used or generated in the project
├── models             # the final conventional and enzyme-constrained GEMs
                       # under different growth conditions,
                       # available in both MATLAB and SBML formats
```

### Installation and Requirements
Ensure that you have [RAVEN](https://github.com/SysBioChalmers/RAVEN), [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html), and the [GUROBI](https://www.gurobi.com/) optimizer installed on your device.

### Citation
```
@article{https://doi.org/10.1111/nph.70528,
author = {Soleymani, Fayaz and Correa, Sandra Marcela and Arend, Marius and Forghanisardaghi, Niayesh and Treves, Haim and Razaghi-Moghadam, Zahra and Nikoloski, Zoran},
title = {Constraint-based metabolic modeling reveals metabolic properties underpinning the unprecedented growth of Chlorella ohadii},
journal = {New Phytologist},
volume = {n/a},
number = {n/a},
pages = {},
keywords = {Chlorella ohadii, de novo model reconstruction, gene targets, genome-scale metabolic model, growth improvement, metabolic model comparison},
doi = {https://doi.org/10.1111/nph.70528},
url = {https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.70528},
eprint = {https://nph.onlinelibrary.wiley.com/doi/pdf/10.1111/nph.70528},
```


