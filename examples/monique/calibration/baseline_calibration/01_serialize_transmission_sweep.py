###########################
# TODO: TO BE CONVERTED   #
###########################

# import os
# import pandas as pd
# import numpy as np
# from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
# from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
# from simtools.SetupParser import SetupParser
# from simtools.ModBuilder import ModBuilder, ModFn
# from malaria.reports.MalariaReport import add_filtered_report, add_summary_report
# import sys
# sys.path.append('../../')
# from simulation.helpers_sim_setup import update_basic_params, habitat_scales, load_master_csv, set_up_hfca
# from simulation.helpers_add_interventions import add_all_interventions
# from simulation.load_paths import load_box_paths
#
# data_path, project_path = load_box_paths(country_name='Example')
#
# expname = 'PfPR_calibration_burnin_sweep_example'
# population_size = 6000  # if this population size is different that what was used for seasonality calibration, need to run '00_rescale_demog_vector_files.py'
# num_seeds = 5
# years = 30
# serialize = True
# pull_from_serialization = False
# itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
#
# # set file describing which intervention inputs are used for this scenario
# scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_for_calibration.csv')
# scen_df = pd.read_csv(scenario_fname)
# # row matching the burnin
# scen_index = scen_df[scen_df['ScenarioName'] == 'transmission_calibration_burnin'].index[0]
#
#
# #############################################################################################################
# # Start from a base MALARIA_SIM config builder, update basic configs, and set up serialization
# cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
# update_basic_params(cb, project_path)
#
# if serialize:
#     cb.update_params( {
#         'Simulation_Duration': years * 365 + 1,
#         'Serialization_Time_Steps': [365*(years-5), 365*years],
#         'Serialized_Population_Reading_Type': 'NONE',
#         'Serialized_Population_Writing_Type': 'TIMESTEP',
#          'Serialization_Mask_Node_Write': 0,  # 0 corresponds to the previous version default: the same larval habitat parameters will be used in the burnin and pickup (from the burnin config)
#          'Serialization_Precision': 'REDUCED'
#         })
# else:
#     cb.update_params({
#         'Serialized_Population_Reading_Type': 'NONE',
#         'Simulation_Duration': years * 365 + 1,
#         'Serialized_Population_Writing_Type': 'NONE'
#     })
#
#
# #############################################################################################################
# # INTERVENTIONS
# # health-seeking
# if (not pd.isna(scen_df.at[scen_index, 'CM_filename'])) and (not (scen_df.at[scen_index, 'CM_filename'] == '')):
#     hs_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'CM_filename']))
# else:
#     hs_df = pd.DataFrame()
# # NMFs
# if (not pd.isna(scen_df.at[scen_index, 'NMF_filename'])) and (not (scen_df.at[scen_index, 'NMF_filename'] == '')):
#     nmf_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'NMF_filename']))
# else:
#     nmf_df = pd.DataFrame()
# # ITNs
# if (not pd.isna(scen_df.at[scen_index, 'ITN_filename'])) and (not (scen_df.at[scen_index, 'ITN_filename'] == '')):
#     itn_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ITN_filename']))
# else:
#     itn_df = pd.DataFrame()
# if (not pd.isna(scen_df.at[scen_index, 'ANC_ITN_filename'])) and (not (scen_df.at[scen_index, 'ANC_ITN_filename'] == '')):
#     itn_anc_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'ANC_ITN_filename']))
# else:
#     itn_anc_df = pd.DataFrame()
# if (not pd.isna(scen_df.at[scen_index, 'EPI_ITN_filename'])) and (not (scen_df.at[scen_index, 'EPI_ITN_filename'] == '')):
#     itn_epi_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'EPI_ITN_filename']))
# else:
#     itn_epi_df = pd.DataFrame()
# itn_decay_params = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'itn_discard_decay_params.csv'))
# itn_use_seasonality = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'ITN_use_seasonality.csv'))
# # IRS
# if (not pd.isna(scen_df.at[scen_index, 'IRS_filename'])) and (not (scen_df.at[scen_index, 'IRS_filename'] == '')):
#     irs_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'IRS_filename']))
# else:
#     irs_df = pd.DataFrame()
# # SMC
# if (not pd.isna(scen_df.at[scen_index, 'SMC_filename'])) and (not (scen_df.at[scen_index, 'SMC_filename'] == '')):
#     smc_df = pd.read_csv(
#         os.path.join(project_path, 'simulation_inputs', '%s.csv' % scen_df.at[scen_index, 'SMC_filename']))
# else:
#     smc_df = pd.DataFrame()
#
#
#
# # CUSTOM REPORTS
# add_filtered_report(cb, start=(years-3)*365, end=years*365)
# add_summary_report(cb, start=(years-1)*365, age_bins=[0.25, 5, 15, 125], interval=30,
#                    description='Monthly', parasitemia_bins=[10, 50, 1e9])
#
# df = load_master_csv(project_path=project_path)
# rel_abundance_df = habitat_scales(project_path=project_path)
# lhdf = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'larval_habitats', 'monthly_habitats_1_%i.csv' % population_size))
#
#
# sweepspace = {
#     # rep_admin: np.logspace(-1.6, 2, 25) for rep_admin in lhdf['archetype'].values
#     rep_admin: [np.logspace(-1.6, 2, 25)[ii] for ii in [0, 4, 8, 12, 16, 20]] for rep_admin in lhdf['archetype'].values   # make sure this matches the values used in the next step
# }
#
# # BUILDER
# builder = ModBuilder.from_list([[ModFn(set_up_hfca, hfca=my_rep_admin,
#                                        archetype_hfca=df.at[my_rep_admin, 'seasonality_archetype'],
#                                        pull_from_serialization=pull_from_serialization,
#                                        ser_date=0,
#                                        hdf=rel_abundance_df,
#                                        lhdf=lhdf,
#                                        population_size=population_size,
#                                        hab_multiplier=hab_scale,
#                                        run_number=x),
#                                  ModFn(add_all_interventions, hfca=my_rep_admin, hs_df=hs_df, nmf_df=nmf_df, itn_df=itn_df,
#                                        itn_anc_df=itn_anc_df, itn_anc_adult_birthday_years=itn_anc_adult_birthday_years,
#                                        itn_epi_df=itn_epi_df, itn_decay_params=itn_decay_params, itn_use_seasonality=itn_use_seasonality,
#                                        irs_df=irs_df, smc_df=smc_df),
#                                  ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
#                                  ModFn(DTKConfigBuilder.set_param, 'Habitat_Multiplier', hab_scale),
#                                  ]
#                                 for my_rep_admin in lhdf['archetype'].values
#                                 for x in range(num_seeds)
#                                 for hab_scale in sweepspace[my_rep_admin]
#                                 ])
#
#
# run_sim_args = {
#     'exp_name': expname,
#     'config_builder': cb,
#     'exp_builder': builder
# }
#
#
# if __name__ == "__main__":
#
#     SetupParser.default_block = 'HPC'
#     SetupParser.init()
#     exp_manager = ExperimentManagerFactory.init()
#     exp_manager.run_simulations(**run_sim_args)
