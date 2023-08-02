# First, run script for a burn-in (burnin_for_seasonalityCalib), set to have a roughly similar transmission intensity as rescaled reference data to establish population immunity
# Second, update the burnin_ids for each archetype (manually enter archetype names) to point to the burn-id, update the rep_admin name, and run the main calibration in two stages:
#   1) run a first round with a longer jump step to get to the general correct area of parameter space (set later_round = False)
#   2) run a second round with shorter jumps to narrow in on the most appropriate parameters (set later_round = True)

import copy
import math
import os
import numpy as np
import pandas as pd
import sys
sys.path.append('../../')
from calibtool.CalibManager import CalibManager
from calibtool.algorithms.OptimTool import OptimTool
from calibtool.plotters.SiteDataPlotter import SiteDataPlotter
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
from simtools.Utilities.COMPSUtilities import COMPS_login
from dtk.vector.species import set_larval_habitat
from malaria.reports.MalariaReport import add_filtered_report
from dtk.interventions.outbreakindividual import recurring_outbreak
from simulation.calibration.helpers_seasonality_calibration import get_spline_values4, \
    get_spline_values4_constantMaxHab, update_starting_spline_values, update_starting_spline_values_constantMaxHab
from simulation.calibration.SeasonalityCalibSite import SeasonalityCalibSite
from simulation.load_paths import load_box_paths
from simulation.helpers_sim_setup import update_basic_params, set_input_files
from simulation.helpers_add_interventions import add_all_interventions

data_path, project_path = load_box_paths(country_name='Example')

itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
simulation_pop = 10000  # from demographics file
vector_human_scalar = simulation_pop / 1000
max_incidence = 100

# specify name of archetype's representative admin
rep_admin = 'AA'
# burn-in id should be listed for each of the archetypes
burnin_ids = {
    'AA': '7c9195e8-3d0a-ee11-aa07-b88303911bc1'
}
# set whether this is the first or second set of calibrations for this archetype
round_number = 1  # <-- user needs to specify which round this is to pick up from prior round. Starts at round_number=1
if round_number > 1:
    later_round = True  # if later_round is True, we use the previous best estimate as the starting point
    if round_number >= 2:
        constant_max_LH = True  # if constant_max_LH is True, we fix maxHab at the best estimate and only adjust other params
    else:
        constant_max_LH = False
else:
    later_round = False
    constant_max_LH = False


pull_from_serialization = True  # should be set to True unless run burn-in within each calibration simulation (warning: that takes a long time)
if pull_from_serialization:
    throwaway = 3  # duration of simulation in years, -1.
    burnin_years = 20  # year when serialization occurred in burnin
    burnin_max_habitat_value = 10 + np.log10(vector_human_scalar)  # make sure this value matches that used in the burnin
    burnin_id = burnin_ids[rep_admin]
else:
    throwaway = 29  # duration of simulation in years, -1.


# set file describing which intervention inputs are used for this scenario
scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_for_calibration.csv')
scen_df = pd.read_csv(scenario_fname)
# row matching this main calibration (not the burnin)
scen_index = scen_df[scen_df['ScenarioName'] == 'seasonality_calibration_noMassITN'].index[0]
expname = '%s_%s_round%i_maxInc%i' % (scen_df.at[scen_index, 'ScenarioName'], rep_admin, round_number,max_incidence)

# set calibration algorithm parameters (some are different between earlier and later round since we make smaller jumps to match target)
r = 0.3  # radius - determines how close the parameter values of each sample within an iteration is to the center value
r_second = 0.005  # 0.1
r_fourth = 0.001
center_move_scale = 1/15
center_move_scale_second = 1/25  # 1/20
center_move_scale_fourth = 1/30
sim_runs_per_param_set = 2  # outside of testing, generally 5
sim_runs_per_param_set_second = 5  # 10
sim_runs_per_param_set_fourth = 10  # 30
max_iterations = 2  # outside of testing, generally 10
max_iterations_second = 20
max_iterations_fourth = 20
samples_per_iteration = 8  # outside of testing, generally 40  # must be at least 8

sites = [SeasonalityCalibSite(hfca=rep_admin, throwaway=throwaway, project_path=project_path)]
plotters = [LikelihoodPlotter(combine_sites=True), SiteDataPlotter(num_to_plot=0, combine_sites=True)]


#############################################################################################################
# Start from a base MALARIA_SIM config builder and update basic configs and serialized file read-in
cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
update_basic_params(cb, project_path)


cb.update_params({
    'Simulation_Duration': (throwaway + 1) * 365,
    'Serialized_Population_Writing_Type': 'NONE',
})

if pull_from_serialization:
    COMPS_login('https://comps.idmod.org')
    print("building from pickup")
    expt = retrieve_experiment(burnin_id)
    # creating data with all the simulation tags
    ser_df = pd.DataFrame([x.tags for x in expt.simulations])
    # getting paths for all the sims
    ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])
    if 'Admin_Name' in ser_df.columns.values:
        ser_df = ser_df[ser_df['MaxHab'] == burnin_max_habitat_value]
        ser_df = ser_df.set_index('Admin_Name')
        ser_path = ser_df.at[rep_admin, 'outpath']
    else:
        ser_path = ser_df['outpath'].values[0]
    cb.update_params({
        'Serialized_Population_Reading_Type': 'READ',
        'Serialized_Population_Path': os.path.join(ser_path, 'output'),
        'Serialized_Population_Filenames': ['state-%05d.dtk' % (burnin_years*365)],
        'Enable_Random_Generator_From_Serialized_Population': 0,
        'Serialization_Mask_Node_Read': 16  # set to 16 so that it will run simulations using larval habitat values from the config file instead of from the burnin
        # 0 corresponds to the previous version default: the same larval habitat parameters will be used in the burnin and pickup (from the burnin config)
    })
else:
    cb.update_params({
        'Serialized_Population_Reading_Type': 'NONE',
    })

#############################################################################################################
# set up the parameters that will be used to calibrate vector larval habitat seasonality

# helper files to specify relative abundance of vector species + initialization of monthly habitats
if constant_max_LH:
    spline, fractions = get_spline_values4_constantMaxHab(rep_admin, project_path)  # for starting guess on monthly multipliers, use rainfall proportional to maximum rainfall month
else:
    spline, fractions = get_spline_values4(rep_admin, project_path)  # for starting guess on monthly multipliers, use rainfall proportional to maximum rainfall month

# if the second round, read in the best values from the previous round and use those as the starting point
if later_round:
    round1_fname = os.path.join('%s_%s_round%i_maxInc%i' % (scen_df.at[scen_index, 'ScenarioName'], rep_admin, (round_number-1), max_incidence), '_plots', 'LL_all.csv')

    round1_df = pd.read_csv(round1_fname)
    # use row with best match to reference dataset
    round1_best_df = round1_df[round1_df.total == round1_df.total.max()].reset_index()
    if constant_max_LH:
        spline = update_starting_spline_values_constantMaxHab(spline, round1_best_df)
    else:
        spline = update_starting_spline_values(spline, round1_best_df)
    if round_number < 4:
        r = r_second
        center_move_scale = center_move_scale_second
        sim_runs_per_param_set = sim_runs_per_param_set_second
        max_iterations = max_iterations_second
    else:
        r = r_fourth
        center_move_scale = center_move_scale_fourth
        sim_runs_per_param_set = sim_runs_per_param_set_fourth
        max_iterations = max_iterations_fourth

params = []
# add months to list of params being sampled
for i, row in spline.iterrows():
    d = row.to_dict()
    params.append(d)

# overall format of monthly habitat to be specified
ls_hab_ref = {'Capacity_Distribution_Number_Of_Years': 1,
              'Capacity_Distribution_Over_Time': {
                  'Times': [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
                  'Values': []
              },
              'Max_Larval_Capacity': 800000000}


# defines how parameters under sampling get consumed by EMOD
def map_sample_to_model_input(cb, sample):

    # set up some variables to hold the values we'll need to update the monthly habitat (hab)
    tags = {}
    sample = copy.deepcopy(sample)

    month_vals = []
    # pull out the variable values at the current step of calibration
    for ii in range(1, 13):
        val_name = 'MonthVal%d' % ii
        cur_value = sample.pop(val_name)
        month_vals.append(cur_value)
        tags.update({val_name: cur_value})
    my_spline = month_vals

    max_habitat_name = 'MaxHab'
    maxvalue = sample.pop(max_habitat_name)
    tags.update({max_habitat_name: maxvalue})
    # constant habitat is currently set as a fixed parameter, not fitted
    const_name = 'Constant'
    const = 3 + np.log10(vector_human_scalar)
    tags.update({const_name: const})

    for (s, sp) in zip(fractions, ['arabiensis', 'funestus', 'gambiae']):
        hab = copy.copy(ls_hab_ref)
        hab['Capacity_Distribution_Over_Time']['Values'] = list(my_spline)
        hab['Max_Larval_Capacity'] = pow(10, maxvalue)*s
        # this function updates EMOD parameters to what is requested based on the calibration parameter sampling
        set_larval_habitat(cb, {sp: {'LINEAR_SPLINE': hab,
                                     'CONSTANT': pow(10, const)*s}})

    # tags.update({'Pop_Scale' : 1})

    return tags


#############################################################################################################
# specify other EMOD setup: climate, demographics, interventions, output reports
set_input_files(cb, rep_admin, rep_admin, population_size=simulation_pop)  # climate and demographics files
recurring_outbreak(cb, start_day=73, outbreak_fraction=0.005, tsteps_btwn=73)

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
add_filtered_report(cb, start=(throwaway-2)*365,
                    end=(throwaway+1)*365)

#############################################################################################################
# calibtool setup
optimtool = OptimTool(params,
                      mu_r=r,           # <-- radius for numerical derivatve.  CAREFUL not to go too small with integer parameters
                      sigma_r=r/10.,  # <-- stdev of radius
                      center_repeats=1,  # <-- Number of times to replicate the center (current guess).  Nice to compare intrinsic to extrinsic noise
                      samples_per_iteration=samples_per_iteration, #  <-- Samples per iteration, includes center repeats.  Actual number of sims run is this number times number of sites.
                      center_move_scale=center_move_scale  # controls
)

calib_manager = CalibManager(name=expname,
                             config_builder=cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=sim_runs_per_param_set,  # <-- Replicates
                             max_iterations=max_iterations,   # <-- Iterations
                             plotters=plotters)

run_calib_args = {
    "calib_manager":calib_manager
}

if __name__ == "__main__":
    # Which simtools.ini block to use for this calibration
    SetupParser.default_block = 'HPC'
    if not SetupParser.initialized:
        SetupParser.init()
    cm = run_calib_args["calib_manager"]
    cm.run_calibration()
