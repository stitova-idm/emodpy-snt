import os
import pandas as pd
import numpy as np
from emodpy_malaria.malaria_config import set_species_param, add_species
from emod_api.interventions.common import change_individual_property_scheduled
from emodpy_malaria.interventions.outbreak import add_outbreak_individual


def update_basic_params(config, manifest, project_path):
    vector_species = ['arabiensis', 'funestus', 'gambiae']
    add_species(config, manifest, vector_species)
    config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    config.parameters.Age_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'
    config.parameters['logLevel_JsonConfigurable'] = 'ERROR'
    config.parameters['logLevel_VectorHabitat'] = 'ERROR'
    config.parameters['logLevel_StandardEventCoordinator'] = 'ERROR'
    config.parameters['Habitat_Multiplier'] = 1
    config.parameters.Enable_Default_Reporting = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Property_Output = 0
    config.parameters.Enable_Demographics_Risk = 1
    config.parameters.Report_Detection_Threshold_Blood_Smear_Parasites = 0
    config.parameters.Report_Parasite_Smear_Sensitivity = 0.01  # number of microliters of blood examined
    config.parameters.Incubation_Period_Distribution = 'CONSTANT_DISTRIBUTION'
    config.parameters.Incubation_Period_Constant = 3  # parameter shortened from 7 to yield a 12-13 day incubation period

    larval_habitat_types_each = {"arabiensis": [{
        "Habitat_Type": "CONSTANT",
        "Max_Larval_Capacity": 10000000.0
    }, {
        "Habitat_Type": "TEMPORARY_RAINFALL",
        "Max_Larval_Capacity": 7500000000.0
    }], "funestus": [{
        "Habitat_Type": "WATER_VEGETATION",
        "Max_Larval_Capacity": 400000000.0
    }], "gambiae": [{
        "Habitat_Type": "CONSTANT",
        "Max_Larval_Capacity": 10000000.0
    }, {
        "Habitat_Type": "TEMPORARY_RAINFALL",
        "Max_Larval_Capacity": 830000000.0
    }]}

    # read in and set vector bionomics
    vector_bionomics = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'vector_bionomics.csv'))
    for vector in vector_species:
        bionomics_row = vector_bionomics[vector_bionomics['species'] == vector].reset_index()
        # overwrite=True is only needed for parameters that are lists (as you can either overwrite or add to them)
        set_species_param(config, vector, 'Anthropophily', bionomics_row['Anthropophily'][0])
        set_species_param(config, vector, 'Indoor_Feeding_Fraction', bionomics_row['Indoor_Feeding_Fraction'][0])
        set_species_param(config, vector, 'Days_Between_Feeds', bionomics_row['Days_Between_Feeds'][0])
        set_species_param(config, vector, 'Habitats', larval_habitat_types_each[vector], overwrite=True)

    # added to match the old config to the new config
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Base_Air_Temperature = 22
    config.parameters.Enable_Vector_Migration = 0


def habitat_scales(project_path):
    rel_abundance_fname = os.path.join(project_path, 'ento', 'DS_vector_rel_abundance.csv')
    rdf = pd.read_csv(rel_abundance_fname)
    rdf = rdf.rename(columns={'Anopheles_arabiensis': 'arabiensis_scale_factor',
                              'Anopheles_coluzzii_gambiae': 'gambiae_scale_factor',
                              'Anopheles_funestus_subgroup': 'funestus_scale_factor'})
    rdf = rdf.set_index('DS')
    return rdf


def set_input_files(config, hfca, archetype_hfca, population_size=1000):
    config.parameters.Climate_Model = "CLIMATE_BY_DATA"
    config.parameters["DS"] = hfca
    config.parameters["Archetype"] = archetype_hfca
    config.parameters.Demographics_Filenames = [os.path.join('demographics_and_climate',
                                                             'season_calib_deomgraphics.json')],
    config.parameters.Air_Temperature_Filename = os.path.join('demographics_and_climate', archetype_hfca,
                                                              '%s_air_temperature_daily_2016.bin' % archetype_hfca),
    config.parameters.Land_Temperature_Filename = os.path.join('demographics_and_climate', archetype_hfca,
                                                               '%s_air_temperature_daily_2016.bin' % archetype_hfca),
    config.parameters.Rainfall_Filename = os.path.join('demographics_and_climate', archetype_hfca,
                                                       '%s_rainfall_daily_2016.bin' % archetype_hfca),
    config.parameters.Relative_Humidity_Filename = os.path.join('demographics_and_climate', archetype_hfca,
                                                                '%s_relative_humidity_daily_2016.bin' % archetype_hfca)

    return {'DS': hfca}


def set_up_hfca(config, platform, hfca, archetype_hfca=None,
                pull_from_serialization=False, ser_date=50 * 365,
                burnin_id='',
                burnin_fname='', population_size=1000,
                hab_multiplier=-1, run_number=-1, use_arch_burnin=False):
    set_input_files(config, hfca, archetype_hfca, population_size)
    if not archetype_hfca:
        archetype_hfca = hfca

    if pull_from_serialization:
        hab_scale_factor_param_name = 'Habitat_Multiplier'
        if burnin_id:
            print("building from pickup: identifying burnin addresses from %s" % burnin_id)
            ser_df = platform.create_sim_directory_df(burnin_id)
        else:
            print("building from pickup: identifying burnin addresses from %s" % burnin_fname)
            ser_df = pd.read_csv(burnin_fname)
        if use_arch_burnin:
            if hab_multiplier >= 0 and run_number >= 0:
                ser_df[hab_scale_factor_param_name] = ser_df[hab_scale_factor_param_name].apply(
                    lambda x: np.round(x, 5))
                sdf = ser_df[(ser_df[hab_scale_factor_param_name] >= (np.round(hab_multiplier, 5) - 0.00001)) &
                             (ser_df[hab_scale_factor_param_name] <= (np.round(hab_multiplier, 5) + 0.00001)) &
                             (ser_df['Run_Number'] == run_number) &
                             (ser_df['admin_name'] == archetype_hfca)]
                config.parameters[hab_scale_factor_param_name] = hab_multiplier
                ser_path = sdf['outpath'].values[0]
            else:
                ser_path = ser_df['outpath'].values[0]
        elif 'admin_name' in ser_df.columns.values:
            ser_df = ser_df[ser_df['Run_Number'] == run_number]
            ser_df = ser_df.set_index('admin_name')
            ser_path = ser_df.at[hfca, 'outpath']
        else:
            ser_path = ser_df['outpath'].values[0]

        config.parameters.Serialized_Population_Reading_Type = 'READ'
        config.parameters.Serialized_Population_Path = os.path.join(ser_path, 'output')
        config.parameters.Serialized_Population_Filenames = ['state-%05d.dtk' % ser_date]
        config.parameters.Enable_Random_Generator_From_Serialized_Population = 0
        config.parameters.Serialization_Mask_Node_Read = 0
        # 0 corresponds to the previous version default: the same larval habitat parameters will be used
        # in the burnin and pickup (from the burnin config)

    else:
        config.parameters.Serialized_Population_Reading_Type = 'NONE'

    return {'admin_name': hfca}


def add_recurring_outbreak(campaign):
    add_outbreak_individual(campaign, demographic_coverage=0.01, start_day=182, repetitions=-1,
                            timesteps_between_repetitions=365)


def load_spline_and_scale_factors(lhdf, archetype_hfca):
    lhdf = lhdf.set_index('archetype')
    my_spline = [lhdf.at[archetype_hfca, 'MonthVal%d' % x] for x in range(1, 13)]
    maxvalue = lhdf.at[archetype_hfca, 'MaxHab']
    const = lhdf.at[archetype_hfca, 'Constant']
    # pop_scale = lhdf.at[archetype_hfca, 'pop_scale']

    return my_spline, maxvalue, const  # , pop_scale


def load_master_csv(project_path):
    master_csv = os.path.join(project_path, 'admin_pop_archetype.csv')
    df = pd.read_csv(master_csv, encoding='latin')
    df['admin_name'] = df['admin_name'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    df = df.set_index('admin_name')
    return df


def update_smc_access_ips(campaign, hfca, smc_df):
    # done
    df = smc_df[smc_df['admin_name'] == hfca]

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # # Original approach was to set IPs at the beginning of the simulation and then update individuals' IPs at their birthday.
    # # However, this does not work as expected because individuals born before the start of the simulation are missed by the birthday-triggered IP change
    # #  (it appears that it is birth-triggered with a delay, so it doesn't get applied to anyone who was already born when the simulation begins).
    #
    # # set IPs at beginning of simulation (does not use IPs from burnin since SMC coverage may change)
    # change_individual_property(cb, 'SMCAccess', 'Low', target={'agemin': 0, 'agemax': 5}, coverage=1, blackout_flag=False)
    # change_individual_property(cb, 'SMCAccess', 'High', target={'agemin': 0, 'agemax': 5}, coverage=df['high_access_U5'].values[0], blackout_flag=False)
    # change_individual_property(cb, 'SMCAccess', 'Low', target={'agemin': 5, 'agemax': 120}, coverage=1, blackout_flag=False)
    # change_individual_property(cb, 'SMCAccess', 'High', target={'agemin': 5, 'agemax': 120}, coverage=df['high_access_5_10'].values[0], blackout_flag=False)
    #
    # # set how IPs change as individuals are born or age
    # change_individual_property_at_age(cb, 'SMCAccess', 'Low', 1, coverage=1)
    # change_individual_property_at_age(cb, 'SMCAccess', 'High', 2, coverage=df['high_access_U5'].values[0])
    # change_individual_property_at_age(cb, 'SMCAccess', 'Low', 365*5, coverage=1)
    # change_individual_property_at_age(cb, 'SMCAccess', 'High', (365*5+1), coverage=df['high_access_5_10'].values[0])
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # change SMCAccess property of newborns (important for those born after the first SMC round in a year)
    if df.shape[0] > 0:
        # SVET - no blackout_flag available, but result is same as blackout_flag=False
        change_individual_property_scheduled(campaign, coverage=1,
                                             new_ip_key='SMCAccess', new_ip_value="Low",
                                             target_age_min=0, target_age_max=5)
        change_individual_property_scheduled(campaign, coverage=df['high_access_U5'].values[0],
                                             new_ip_key='SMCAccess', new_ip_value="High",
                                             target_age_min=0, target_age_max=5)

        # before the first SMC round in each year, change the SMCAccess IP for the U5 and O5 age groups
        # simdays of the first rounds (change IPs one week before)
        first_round_days = df.loc[df['round'] == 1, 'simday'].values
        change_ip_days = [first_round_days[yy] - 7 for yy in range(len(first_round_days))]

        for rr in change_ip_days:
            change_individual_property_scheduled(campaign, start_day=rr, coverage=1,
                                                 new_ip_key='SMCAccess', new_ip_value="Low",
                                                 target_age_min=0, target_age_max=5)
            change_individual_property_scheduled(campaign, start_day=rr, coverage=df['high_access_U5'].values[0],
                                                 new_ip_key='SMCAccess', new_ip_value="High",
                                                 target_age_min=0, target_age_max=5)
            change_individual_property_scheduled(campaign, start_day=rr, coverage=1,
                                                 new_ip_key='SMCAccess', new_ip_value="Low",
                                                 target_age_min=5, target_age_max=120)
            change_individual_property_scheduled(campaign, start_day=rr, coverage=df['high_access_5_10'].values[0],
                                                 new_ip_key='SMCAccess', new_ip_value="High",
                                                 target_age_min=5, target_age_max=120)

    return {'admin_name': hfca}


# SVET - I can't find what this is doing or where it's being used..
# but this should be straightforward to implement whenever this needs to be used
# def update_drug_config(cb, list_drug_param, value, list_drug_param2=[], value2=0, list_drug_param3=[], value3=0,
#                        list_drug_param4=[], value4=0):
#     cb.update_params({
#         "Malaria_Drug_Params": update_drugs(list_drug_param=list_drug_param, value=value,
#                                             list_drug_param2=list_drug_param2, value2=value2,
#                                             list_drug_param3=list_drug_param3, value3=value3,
#                                             list_drug_param4=list_drug_param4, value4=value4, )
#     })
#     return ({'new drug value': value})


if __name__ == "__main__":
    pass
