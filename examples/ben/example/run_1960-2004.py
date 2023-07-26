import os
import pandas as pd
import time
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import update_species_param
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from load_paths import load_paths
from hbhi.set_up_general import initialize_cb, setup_ds
from tqdm import tqdm

SetupParser.default_block = 'NUCLUSTER'
iopath = load_paths()
larval_hab_csv = 'simulation_inputs/monthly_habitats.csv'
master_csv = 'simulation_inputs/country_DS.csv'
rel_abund_csv = 'simulation_inputs/DS_vector_rel_abundance.csv'

homepath = os.path.expanduser('~')
user = homepath.split('/')[2]
expname = f'{user}_dtknu_1960-2004'

num_seeds = 1
years = 45
serialize = True
pull_from_serialization = False


# BASIC SETUP
# Filtered report for the last year
cb = initialize_cb(years, serialize, filtered_report=1)
cb.update_params( {
    'x_Temporary_Larval_Habitat': 1, # Package default is 0.2
    'x_Base_Population': 1,
    'x_Birth': 1
})

# Changing species param from defaults
update_species_param(cb, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
update_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
update_species_param(cb, 'funestus', 'Anthropophily', 0.5, overwrite=True)
update_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
update_species_param(cb, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)

# Important DFs
master_df = pd.read_csv(os.path.join(iopath, master_csv))
master_df = master_df.set_index('DS_Name')
lhdf = pd.read_csv(os.path.join(iopath, larval_hab_csv))
rel_abund_df = pd.read_csv(os.path.join(iopath, rel_abund_csv))
rel_abund_df = rel_abund_df.set_index('DS_Name')

samp_df = pd.read_csv(os.path.join(iopath, 'simulation_priors/selected_particles.csv'))
ds_list = master_df.archetype.unique()


# BUILDER
list_of_sims = []
for my_ds in tqdm(ds_list):
    samp_ds = samp_df[samp_df.archetype == my_ds]
    # For each archetype, each Habitat_Multiplier has same reduce_id and seed1 by design of preprocessing
    samp_ds1 = samp_ds[['Habitat_Multiplier', 'reduce_id', 'seed1']].drop_duplicates()
    L = []
    for r, row in samp_ds1.iterrows():
        L = L + [[ModFn(setup_ds,
                        my_ds=my_ds,
                        archetype_ds=master_df.at[my_ds, 'archetype'],
                        pull_from_serialization=pull_from_serialization,
                        burnin_id='',
                        ser_date=0,
                        rel_abund_df=rel_abund_df,
                        lhdf=lhdf,
                        demographic_suffix='_wSMC_risk_wIP',
                        climate_prefix=False,
                        use_arch_burnin=False,
                        use_arch_input=True,
                        hab_multiplier=row['Habitat_Multiplier'],
                        parser_default=SetupParser.default_block),
                  ModFn(DTKConfigBuilder.set_param, 'Habitat_Multiplier', row['Habitat_Multiplier']),
                  ModFn(DTKConfigBuilder.set_param, 'Run_Number', row['seed1']),
                  ModFn(DTKConfigBuilder.set_param, 'Sample_ID', row['reduce_id'])]]
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
