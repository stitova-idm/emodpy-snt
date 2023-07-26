import os
import pandas as pd
import time
from scipy import interpolate
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import update_species_param
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from load_paths import load_paths
from hbhi.set_up_general import initialize_cb, setup_ds
from hbhi.set_up_interventions import InterventionSuite, add_all_interventions, update_smc_access_ips
from hbhi.utils import (add_monthly_parasitemia_rep_by_year, add_nmf_trt, tryread_df,
                        add_annual_parasitemia_rep)
from tqdm import tqdm


SetupParser.default_block = 'NUCLUSTER'
iopath = load_paths()
larval_hab_csv = 'simulation_inputs/monthly_habitats.csv'
master_csv = 'simulation_inputs/country_DS.csv'
rel_abund_csv = 'simulation_inputs/DS_vector_rel_abundance.csv'
samp_csv = os.path.join(iopath, 'simulation_priors/selected_particles.csv')

homepath = os.path.expanduser('~')
user = homepath.split('/')[2]
expname = f'{user}_dtknu_2005-2022'

num_seeds = 3
years = 18  # 2005 to 2022 year
ser_num_seeds = 1
ser_date = 45*365
serialize = True
pull_from_serialization = True
use_arch_burnin = True

burnin_id = "2023_06_13_20_38_23_696396"

homepath = os.path.expanduser('~')
user = homepath.split('/')[2]



# BASIC SETUP
cb = initialize_cb(years, serialize=serialize, filtered_report=years,
                   event_reporter=False)
cb.update_params({
    'x_temporary_Larval_Habitat': 1, # Package default is 0.2
    'x_Base_Population': 1,
    'x_Birth': 1
})
cb.update_params({
    "Report_Event_Recorder": 1,
    "Report_Event_Recorder_Individual_Properties": [],
    "Report_Event_Recorder_Events" : ['Received_NMF_Treatment', 'Received_Severe_Treatment',
                                      'Received_Treatment'],
    "Report_Event_Recorder_Ignore_Events_In_List": 0
})
update_species_param(cb, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
update_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
update_species_param(cb, 'funestus', 'Anthropophily', 0.5, overwrite=True)
update_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
update_species_param(cb, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)


# INTERVENTIONS
# Intervention Suite
int_suite =  InterventionSuite()
# health seeking
int_suite.hs_ds_col = 'DS_Name'
int_suite.hs_duration = 365 # Set to none when want to follow df's duration to the dot
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

# Treatment-seeking and malaria testing prompted by NMF
add_nmf_trt(cb, years, 0)

# Treatment-seeking
hs_df = tryread_df(os.path.join(iopath, 'simulation_inputs',
                                'scenarios', 'cm_2005-2022.csv'))
# ITNs
itn_df = tryread_df(os.path.join(iopath, 'simulation_inputs',
                                'scenarios', 'itn_2005-2022.csv'))
# SMC
smc_df = tryread_df(os.path.join(iopath, 'simulation_inputs',
                                'scenarios', 'smc_2005-2022.csv'))
smc_df['high_access_5_10'] = smc_df['high_access_O5'] # update_smc_access_ips hard-coded the column name...
smc_df['coverage_low_access_U5'] = 0
smc_df['coverage_low_access_O5'] = 0


# CUSTOM REPORTS
add_monthly_parasitemia_rep_by_year(cb, num_year=years, tot_year=years,
                                    sim_start_year=2005,
                                    yr_plusone=True, prefix='Monthly_')
add_annual_parasitemia_rep(cb, num_year=years, tot_year=years,
                           sim_start_year=2005,
                           age_bins = [0.25, 1, 2, 5, 10, 15, 30, 50, 125])

# Important DFs
master_df = pd.read_csv(os.path.join(iopath, master_csv))
master_df = master_df.set_index('DS_Name')
lhdf = pd.read_csv(os.path.join(iopath, larval_hab_csv))
rel_abund_df = pd.read_csv(os.path.join(iopath, rel_abund_csv))
rel_abund_df = rel_abund_df.set_index('DS_Name')

samp_df = pd.read_csv(samp_csv)
ds_list = samp_df.DS_Name.unique()

# BUILDER
def lin_interpolate(years, vals):
    val_interp = interpolate.interp1d(years, vals)
    full_vals = val_interp(range(years[0], years[-1]+1))
    return (full_vals)

list_of_sims = []
for my_ds in tqdm(ds_list):
    samp_ds = samp_df[samp_df.DS_Name == my_ds].copy()
    L = []
    for r, row in samp_ds.iterrows():
        hs_ds = hs_df.copy()[hs_df[int_suite.hs_ds_col] == my_ds]
        full_covs = lin_interpolate([2005, 2013, 2017, 2020, 2022],
                                    [0, row['CM_cov_2013'], row['CM_cov_2017'],
                                     row['CM_cov_2020'], row['CM_cov_2020']])
        hs_ds['U5_coverage'] = full_covs * hs_ds['U5_coverage']
        hs_ds['adult_coverage'] = full_covs * hs_ds['adult_coverage']

        itn_ds = itn_df.copy()[itn_df[int_suite.itn_ds_col] == my_ds]
        itn_ds = itn_ds[itn_ds['year'] <= 2022]
        itn_ds = itn_ds.reset_index()

        # TODO: Make code below more defensive
        itn_ds.loc[0, int_suite.itn_cov_cols] = itn_ds.loc[0, int_suite.itn_cov_cols] * row['ITN_2011']
        itn_ds.loc[1, int_suite.itn_cov_cols] = itn_ds.loc[1, int_suite.itn_cov_cols] * row['ITN_2014']
        itn_ds.loc[2, int_suite.itn_cov_cols] = itn_ds.loc[2, int_suite.itn_cov_cols] * row['ITN_2017']
        itn_ds.loc[3, int_suite.itn_cov_cols] = itn_ds.loc[3, int_suite.itn_cov_cols] * row['ITN_2020']

        L = L + [[ModFn(setup_ds,
                        my_ds=my_ds,
                        archetype_ds=master_df.at[my_ds, 'archetype'],
                        pull_from_serialization=pull_from_serialization,
                        burnin_id=burnin_id,
                        ser_date=ser_date,
                        rel_abund_df=rel_abund_df,
                        lhdf=lhdf,
                        demographic_suffix='_wSMC_risk_wIP',
                        climate_prefix=False,
                        use_arch_burnin=True,
                        use_arch_input=False,
                        hab_multiplier=row['Habitat_Multiplier'],
                        parser_default=SetupParser.default_block,
                        serialize_match_tag=['Habitat_Multiplier'],
                        serialize_match_val=[float(row['Habitat_Multiplier'])]),
                   ModFn(add_all_interventions,
                         int_suite=int_suite,
                         my_ds=my_ds,
                         hs_df=hs_ds,
                         itn_df=itn_ds,
                         smc_df=smc_df,
                         addtl_smc_func=update_smc_access_ips, # Change IP every year
                         rep_start=0,
                         rep_duration=years*365),
                  ModFn(DTKConfigBuilder.set_param, 'DS_Name', my_ds),
                  ModFn(DTKConfigBuilder.set_param, 'Sample_ID', row['id']),
                  ModFn(DTKConfigBuilder.set_param, 'Run_Number', row['seed2']+x)]
                  for x in range(num_seeds)]
    list_of_sims = list_of_sims + L

builder = ModBuilder.from_list(list_of_sims)

run_sim_args = {
    'exp_name': expname,
    'config_builder': cb,
    'exp_builder': builder
}

if __name__ == "__main__":
    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)

    time.sleep(20)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
    expt_id = exp_manager.experiment.exp_id

    print(f"python simulation/analyzer/analyze_2005-2022.py -name 2005-2022 -id {expt_id}")
    os.system(f"python simulation/analyzer/analyze_2005-2022.py -name 2005-2022 -id {expt_id}")

    #print(f"Rscript output_processor/iptp_postprocessing.R past ")
    #os.system(f"Rscript output_processor/iptp_postprocessing.R past")
