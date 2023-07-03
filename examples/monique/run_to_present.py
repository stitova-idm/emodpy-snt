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
# from simulation.load_paths import load_box_paths
# from simulation.helpers_add_interventions import add_all_interventions
# from simulation.helpers_sim_setup import update_basic_params, load_master_csv, habitat_scales, set_up_hfca, get_burnin_exp
#
#
# data_path, project_path = load_box_paths(country_name='Example')
#
# num_seeds = 50
# start_year = 2010
# years = 12  # beginning of 2010 to beginning of 2010+ year
# serialize = True
# pull_from_serialization = True
# num_burnin_seeds = 5  # number of seeds run during previous burnin and habitat multiplier calibration
# ser_date = 30*365
# itn_anc_adult_birthday_years = [20, 22, 24, 26, 28]
# population_size = 6000  # needs to match burnin simulation population size
#
# use_arch_burnin = True
# burnin_id = 'cab2ccbc-f70a-ee11-aa07-b88303911bc1'  # <- on IDMCloud  # generated from serialize_transmission_sweep (1960-2010)
#
# scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_to_present.csv')
# scen_df = pd.read_csv(scenario_fname)
# scen_index = scen_df[scen_df['status'] == 'run'].index[0]
# expname = scen_df.at[scen_index, 'ScenarioName']
#
#
# #############################################################################################################
# # Start from a base MALARIA_SIM config builder, update basic configs, set up serialization pickup, and add interventions
# cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
# update_basic_params(cb, project_path)
# ser_df = get_burnin_exp(burnin_id=burnin_id)
#
# cb.update_params({
#         "Report_Event_Recorder": 1,
#         "Report_Event_Recorder_Events": ['Received_Severe_Treatment'],  #!, 'Bednet_Got_New_One'],
#         "Report_Event_Recorder_Ignore_Events_In_List": 0,
#         "Report_Event_Recorder_Individual_Properties": []
#         })
# if serialize:
#     cb.update_params({
#         'Simulation_Duration': years * 365 + 1,
#         'Serialization_Time_Steps': [365*(years-5), 365*years],
#         'Serialized_Population_Reading_Type': 'NONE',
#         'Serialized_Population_Writing_Type': 'TIMESTEP',
#         'Serialization_Mask_Node_Write': 0,  # 0 corresponds to the previous version default: the same larval habitat parameters will be used in the burnin and pickup (from the burnin config)
#         'Serialization_Precision': 'REDUCED',
#         })
# else:
#     cb.update_params({
#         'Serialized_Population_Reading_Type': 'NONE',
#         'Simulation_Duration': years * 365 + 1,
#         'Serialized_Population_Writing_Type': 'NONE'
#     })
#
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
# # CUSTOM REPORTS
# add_filtered_report(cb, start=0, end=years * 365)
# for year in range(years):
#     add_summary_report(cb, start=365 * year, age_bins=[0.25, 5, 15, 30, 50, 125], interval=30, duration_days=365,
#                        description='Monthly%d' % (year + start_year), parasitemia_bins=[10, 50, 1e9])
#     # add_summary_report(cb, start=365 * year, age_bins=[1, 5, 120], interval=30, duration_days=365,
#     #                    description='Monthly_U1U5_%d' % (year + start_year), parasitemia_bins=[10, 50, 1e9])
#
# # FOR CONFIGURING LARVAL HABTIATS
# df = load_master_csv(project_path=project_path)
# hab_scale_factor_fname = os.path.join(project_path, 'simulation_inputs', 'larval_habitats', 'larval_habitat_multipliers_v1.csv')
# hab_df = pd.read_csv(hab_scale_factor_fname)
# rel_abundance_df = habitat_scales(project_path=project_path)
# lhdf = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'larval_habitats', 'monthly_habitats_1_%i.csv' % population_size))
#
# # BUILDER
# builder = ModBuilder.from_list([[ModFn(set_up_hfca, hfca=my_admin,
#                                        archetype_hfca=df.at[my_admin, 'seasonality_archetype'],
#                                        pull_from_serialization=pull_from_serialization,
#                                        ser_date=ser_date,
#                                        hdf=rel_abundance_df,
#                                        lhdf=lhdf,
#                                        population_size=population_size,
#                                        # get the habitat multiplier that matches this admin and is equal to this seed modulo num_burnin_seeds
#                                        hab_multiplier=(hab_df.loc[np.logical_and.reduce([hab_df[c] == v for c, v in
#                                             zip(['admin_name', 'Run_Number'], [my_admin, (x % num_burnin_seeds)])])].reset_index().at[0, 'Habitat_Multiplier']),
#                                        run_number=(x % num_burnin_seeds),
#                                        use_arch_burnin=use_arch_burnin,
#                                        ser_df=ser_df),
#                                 ModFn(add_all_interventions, seed_index=x+1, hfca=my_admin, hs_df=hs_df, nmf_df=nmf_df,
#                                       itn_df=itn_df,
#                                       itn_anc_df=itn_anc_df, itn_anc_adult_birthday_years=itn_anc_adult_birthday_years,
#                                       itn_epi_df=itn_epi_df,
#                                       itn_decay_params=itn_decay_params,
#                                       itn_use_seasonality=itn_use_seasonality,
#                                       irs_df=irs_df, smc_df=smc_df),
#
#                                  ModFn(DTKConfigBuilder.set_param, 'Run_Number', x),
#                                  ModFn(DTKConfigBuilder.set_param, 'Habitat_Multiplier',
#                                        (hab_df.loc[np.logical_and.reduce([hab_df[c] == v for c, v in
#                                         zip(['admin_name', 'Run_Number'], [my_admin, (x % num_burnin_seeds)])])].reset_index().at[0, 'Habitat_Multiplier'])),
#                                  ]
#                                 # for my_admin in ['Gitega']
#                                 # for my_admin in ['Fota', 'Gitega']
#                                 for my_admin in df.index
#                                 for x in range(num_seeds)
#                                 ])
#
# run_sim_args = {
#     'exp_name': expname,
#     'config_builder': cb,
#     'exp_builder': builder
# }
#
# if __name__ == "__main__":
#
#     SetupParser.default_block = 'HPC'
#     SetupParser.init()
#     exp_manager = ExperimentManagerFactory.init()
#     exp_manager.run_simulations(**run_sim_args)
#
