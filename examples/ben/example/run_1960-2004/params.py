import os
import manifest
import pandas as pd
from pathlib import Path
from snt.hbhi.load_paths import load_paths

test_run = False

iopath = load_paths(user_path=manifest.BEN_DIR)
larval_hab_csv = 'simulation_inputs/monthly_habitats.csv'
master_csv = 'simulation_inputs/country_DS.csv'
rel_abund_csv = 'simulation_inputs/DS_vector_rel_abundance.csv'
samp_csv = os.path.join(iopath, 'simulation_priors/selected_particles.csv')

homepath = os.path.expanduser('~')
user = Path(homepath).name
expname = f'{user}_dtknu_1960-2004'

num_seeds = 1
years = 45
serialize = True
pull_from_serialization = False

# Add new
filtered_report = 1  # represent  num_year
yr_plusone = True
event_reporter = False

# ZDU add
demographic_suffix = '_wSMC_risk_wIP'
climate_suffix = ''
climate_prefix = False
use_arch_burnin = True
use_arch_input = True  # ZDU: test False

rep_start = 0
rep_duration = 10000

# load dfs
master_df = pd.read_csv(os.path.join(iopath, master_csv))
master_df = master_df.set_index('DS_Name')
ds_list = master_df.archetype.unique()

samp_df = pd.read_csv(samp_csv)
samp_df = samp_df[samp_df.DS_Name.isin(['Anie', 'Doufelgou'])]
