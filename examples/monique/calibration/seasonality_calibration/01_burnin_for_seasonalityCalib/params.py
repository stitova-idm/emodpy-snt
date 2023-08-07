import os
import pandas as pd
import numpy as np
import manifest

rep_admin = 'AA'
simulation_pop = 10000  # from demographics file
vector_human_scalar = simulation_pop / 1000

max_habitat_value = 10 + np.log10(vector_human_scalar)
month_scalar = 0.01
const_habitat = 3 + np.log10(vector_human_scalar)
num_seeds = 5
years = 20
serialize = True
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]

# ZDU
demographics_file = os.path.join('demographics_and_climate', '_entire_country',
                                 f'demographics_each_admin_{simulation_pop}.json')

# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(manifest.project_path, 'simulation_inputs', '_intervention_file_references',
                              'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching this main calibration (not the burnin)
scen_index = scen_df[scen_df['ScenarioName'] == 'seasonality_calibration_burnin'].index[0]
expname = '%s_%s_maxHab%d_monthScalar%d' % (
    scen_df.at[scen_index, 'ScenarioName'], rep_admin, max_habitat_value * 1000, month_scalar * 1000)
