import os
import pandas as pd
from pathlib import Path
from snt.hbhi.load_paths import load_paths
from snt.hbhi.set_up_interventions import InterventionSuite
import manifest

test_run = False

iopath = load_paths(user_path=manifest.BEN_DIR)
larval_hab_csv = 'simulation_inputs/monthly_habitats.csv'
master_csv = 'simulation_inputs/country_DS.csv'
rel_abund_csv = 'simulation_inputs/DS_vector_rel_abundance.csv'
samp_csv = os.path.join(iopath, 'simulation_priors/selected_particles.csv')

homepath = os.path.expanduser('~')
user = Path(homepath).name
expname = f'{user}_dtknu_2005-2022'

num_seeds = 1  #
years = 18  # 2005 to 2022 year
ser_num_seeds = 1
ser_date = 45 * 365
serialize = True
pull_from_serialization = True
use_arch_burnin = True

# Add new
demographic_suffix = '_wSMC_risk_wIP'
climate_suffix = ''
climate_prefix = False
use_arch_input = False
yr_plusone = True
event_reporter = False
add_event_report = True
filtered_report = 1  # represent  num_year

burnin_id = "2023_07_27_13_50_24_524201"

# INTERVENTIONS
# Intervention Suite
int_suite = InterventionSuite()
# health seeking
int_suite.hs_ds_col = 'DS_Name'
int_suite.hs_duration = 365  # Set to none when want to follow df's duration to the dot
# itn
int_suite.itn_ds_col = 'DS_Name'
int_suite.itn_discard_distribution = 'weibull'
int_suite.itn_cov_cols = ['U05', 'U10', 'U20', 'A20']
int_suite.itn_cov_age_bin = [0, 5, 10, 20]
int_suite.itn_retention_in_yr = 1.69
int_suite.itn_seasonal_months = [0, 91, 182, 274]
int_suite.itn_seasonal_values = [0.739, 0.501, 0.682, 1]
# smc
int_suite.smc_adherence = True
int_suite.smc_coverage_col = ['coverage_high_access_U5', 'coverage_low_access_U5',
                              'coverage_high_access_O5', 'coverage_low_access_O5']
int_suite.smc_access = ['High', 'Low', 'High', 'Low']
int_suite.smc_agemins = [0.25, 0.25, 5, 5]
int_suite.smc_agemax_type = ['fixed', 'fixed', 'fixed', 'fixed']
int_suite.smc_agemaxs = [5, 5, 10, 10]
int_suite.smc_leakage = False
# irs unchanged

# load dfs
master_df = pd.read_csv(os.path.join(iopath, master_csv))
master_df = master_df.set_index('DS_Name')

samp_df = pd.read_csv(samp_csv)
samp_df = samp_df[samp_df.DS_Name.isin(['Anie', 'Doufelgou'])]
ds_list = samp_df.archetype.unique()
