# Run Snakemake Examples #
1) Create virtual environment in emodpy-snt
```
PATH_TO_EMODPY_SNT\emodpy-snt\python -m venv <location_and_name_of_venv>
```

2) Activate environment
```
<location_and_name_of_venv>\Scripts\activate.bat
```

3) pip
```
venv> pip install emodpy_snt --index-url=https://packages.idmod.org/api/pypi/pypi-production/simple

```

4) Create profile for snakemake (see snakemake help)

On Windows:
In C:\ProgramData\snakemake\snakemake or C:\Users\USER_NAME\AppData\Local\snakemake\snakemake create a directory ```default``` and a file ```default.yaml``` (or ```config.yaml```, depending on machine's ask)
with the content: 
```
jobs: 1
```

5) Delete existing output files
```
venv> snakemake --profile=default --cores=2 clean
```

6) Run all snakemake rules in _snakefile_
```
venv> snakemake --profile=default --cores=2
```

To run only one rule from the set of rules in _snakefile_
```
venv> snakemake --profile=default --cores=2 snakemake_rule
```

