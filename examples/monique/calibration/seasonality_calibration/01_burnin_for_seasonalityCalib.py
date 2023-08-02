# First part of calibration process. Should be set to have a roughly similar transmission intensity as rescaled reference data to establish population immunity

import sys
import os
sys.path.append('../../')
import copy
import pandas as pd
import numpy as np
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_larval_habitat
from dtk.interventions.outbreakindividual import recurring_outbreak
from simtools.SetupParser import SetupParser
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from malaria.reports.MalariaReport import add_filtered_report
from simulation.calibration.helpers_seasonality_calibration import get_burnin_spline_values
from simulation.helpers_sim_setup import update_basic_params, set_input_files
from simulation.helpers_add_interventions import add_all_interventions
from simulation.load_paths import load_box_paths

data_path, project_path = load_box_paths(country_name='Example')

# specify archetype representative name and the maximum larval habitat for the burnin.
rep_admin = 'AA'
simulation_pop = 10000  # from demographics file
vector_human_scalar = simulation_pop / 1000

max_habitat_value = 10 + np.log10(vector_human_scalar)
month_scalar = 0.01
const_habitat = 3 + np.log10(vector_human_scalar)
num_seeds = 1
years = 20
serialize = True
itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]


# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching this main calibration (not the burnin)
scen_index = scen_df[scen_df['ScenarioName'] == 'seasonality_calibration_burnin'].index[0]
expname = '%s_%s_maxHab%d_monthScalar%d' % (scen_df.at[scen_index, 'ScenarioName'], rep_admin, max_habitat_value*1000, month_scalar*1000)

#############################################################################################################
# Start from a base MALARIA_SIM config builder, update basic configs, and set up serialization
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
update_basic_params(cb, project_path)

if serialize :
    cb.update_params( {
        'Simulation_Duration': years * 365 + 1,
        'Serialization_Time_Steps' : [365*(years-5), 365*years],
        'Serialized_Population_Reading_Type': 'NONE',
        'Serialized_Population_Writing_Type': 'TIMESTEP',
         'Serialization_Mask_Node_Write': 0,  # 0 corresponds to the previous version default: the same larval habitat parameters will be used in the burnin and pickup (from the burnin config)
         'Serialization_Precision': 'REDUCED'
        })
else:
    cb.update_params({
        'Serialized_Population_Reading_Type': 'NONE',
        'Simulation_Duration': years * 365 + 1,
        'Serialized_Population_Writing_Type': 'NONE'
    })


#############################################################################################################
# specify vectors and monthly habitat
month_vals, fractions = get_burnin_spline_values(rep_admin, month_scalar, project_path)
ls_hab_ref = {'Capacity_Distribution_Number_Of_Years': 1,
              'Capacity_Distribution_Over_Time': {
                   'Times' : [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
                   'Values' : month_vals
               },
              'Max_Larval_Capacity' : pow(10, max_habitat_value)}

for (s, sp) in zip(fractions, ['arabiensis', 'funestus', 'gambiae']):
    hab = copy.copy(ls_hab_ref)
    hab['Max_Larval_Capacity'] = pow(10, max_habitat_value)*s
    set_larval_habitat(cb, {sp: {'LINEAR_SPLINE': hab,
                                 'CONSTANT': pow(10, const_habitat)*s}})


#############################################################################################################
# other EMOD setup: climate, demographics, interventions, output reports
set_input_files(cb, rep_admin, rep_admin, population_size=simulation_pop)  # climate and demographics files
recurring_outbreak(cb, start_day=182, outbreak_fraction=0.01, tsteps_btwn=365)

# INTERVENTIONS
# health-seeking
if (not pd.isna(scen_df.at[scen_index, 'CM_filename'])) and (not (scen_df.at[scen_index, 'CM_filename'] == '')):
    hs_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'CM_filename']))
else:
    hs_df = pd.DataFrame()
# NMFs
if (not pd.isna(scen_df.at[scen_index, 'NMF_filename'])) and (not (scen_df.at[scen_index, 'NMF_filename'] == '')):
    nmf_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'NMF_filename']))
else:
    nmf_df = pd.DataFrame()
# ITNs
if (not pd.isna(scen_df.at[scen_index, 'ITN_filename'])) and (not (scen_df.at[scen_index, 'ITN_filename'] == '')):
    itn_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ITN_filename']))
else:
    itn_df = pd.DataFrame()
if (not pd.isna(scen_df.at[scen_index, 'ANC_ITN_filename'])) and (not (scen_df.at[scen_index, 'ANC_ITN_filename'] == '')):
    itn_anc_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ANC_ITN_filename']))
else:
    itn_anc_df = pd.DataFrame()
if (not pd.isna(scen_df.at[scen_index, 'EPI_ITN_filename'])) and (not (scen_df.at[scen_index, 'EPI_ITN_filename'] == '')):
    itn_epi_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'EPI_ITN_filename']))
else:
    itn_epi_df = pd.DataFrame()
itn_decay_params = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'itn_discard_decay_params.csv'))
itn_use_seasonality = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'ITN_use_seasonality.csv'))
# IRS
if (not pd.isna(scen_df.at[scen_index, 'IRS_filename'])) and (not (scen_df.at[scen_index, 'IRS_filename'] == '')):
    irs_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'IRS_filename']))
else:
    irs_df = pd.DataFrame()
# SMC
if (not pd.isna(scen_df.at[scen_index, 'SMC_filename'])) and (not (scen_df.at[scen_index, 'SMC_filename'] == '')):
    smc_df = pd.read_csv(
        os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'SMC_filename']))
else:
    smc_df = pd.DataFrame()

add_all_interventions(cb, hfca=rep_admin, hs_df=hs_df, nmf_df=nmf_df, itn_df=itn_df,
                      itn_anc_df=itn_anc_df, itn_anc_adult_birthday_years=itn_anc_adult_birthday_years,
                      itn_epi_df=itn_epi_df, itn_decay_params=itn_decay_params, itn_use_seasonality=itn_use_seasonality,
                      irs_df=irs_df, smc_df=smc_df)


# request output files
add_filtered_report(cb, start=(years-5)*365, end=years*365)
# add_filtered_report(cb, start=0, end=years*365)


# BUILDER
builder = ModBuilder.from_list([[ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
                                 ModFn(DTKConfigBuilder.set_param, 'MaxHab', max_habitat_value),
                                 ModFn(DTKConfigBuilder.set_param, 'ConstHab', const_habitat),
                                 ModFn(DTKConfigBuilder.set_param, 'Admin_Name', my_rep_admin),
                                 ]
                                for my_rep_admin in [rep_admin]
                                for x in range(num_seeds)
                                ])

run_sim_args = {
    'exp_name': expname,
    'config_builder': cb,
    'exp_builder': builder
}


if __name__ == "__main__":

    SetupParser.default_block = 'HPC'
    SetupParser.init()
    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
