# snt-dtknu

This repository contains a working example of SNT modeling using dtk-tools. It is modified from hbhi-togo, and contains all the input files (CSVs) necessary to run the example. The code is meant to be used in an `NUCLUSTER` environment.

## Explanation

In this example, we have already completed the calibration process. For each district (DS), we have selected the best 5 "particles", or the parameter sets, that fit to the empirical data. So the idea here is based on these selected particles, we go back to the beginning to run the "burnin" (1960-2004), and then pick up to run from 2005 to 2022.

Each district is assigned to an archetype. For the burnin, we only run simulation for the archetypes. In this example there's only one archetype `Doufelgou`, but in actual simulations there will be more archetypes. Then during the pickup, we use district-specific climates, and pull out the serialized burnin simulation that contains the matching archetype and habitat multiplier and run the simulation from there. Each district has 5 selected particles, and for each particles we run 3 realizations (seeds).


## Layout

    .
    ├── IO                      # Simulations input and output
    │   ├── simulation_inputs   # Simulation inputs (climate files, csvs for interventions etc.)
    │   └── simulation_priors   # Contains parameters
    ├── simulation              # Codes related to running the simulations
    │   ├── analyzer            # Analyzing scripts and the collection of pre-written analyzers
    │   └── example             # Codes for running a burnin and a pickup
    ├── simtools.ini            # Do change the paths inside simtools.ini
    ├── load_paths.py           # Do change the paths here pointing to the IO folder
    └── README.md

## Running the code example

You should run the code example from the base directory:
```
python simulation/example/run_1960-2004.py
python simulation/example/run_2005-2022.py
```

