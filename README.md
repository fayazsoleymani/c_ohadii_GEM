## Genome-scale metabolic reconstruction (GEM) of *Chlorella ohadii*

Here, we present a fully automated platform for the *de novo* generation of genome-scale metabolic models (GEMs). We deployed this platform to reconstruct the GEM and the enzyme-constrained GEM for *C. ohadii*.

### Contents
Here, the project's directory structure is illustrated: 
```
c_ohadii_GEM
├── code
│   ├── functions      # contains the implementation of the main functions
│   └── examples.m     # provides examples of how to use the functions
├── data               # the data, used or generated in the project
├── models             # the final conventional and enzyme-constrained GEMs
                       # under different growth conditions,
                       # available in both MATLAB and SBML formats
```

### Installation and Requirements
Ensure that you have [RAVEN](https://github.com/SysBioChalmers/RAVEN), [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/index.html), and the [GUROBI](https://www.gurobi.com/) optimizer installed on your device.

### Citation
TODO


