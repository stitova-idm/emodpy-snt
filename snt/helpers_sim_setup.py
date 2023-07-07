import os
import pandas as pd
import numpy as np
from emodpy_malaria.interventions.outbreak import add_outbreak_individual
from emodpy_malaria.malaria_config import configure_linear_spline, set_species_param, add_species
from simtools.Utilities.Experiments import retrieve_experiment, retrieve_simulation
from emod_api.interventions.common import change_individual_property_scheduled
from dtk.utils.reports.CustomReport import BaseReport
from COMPS.Data import Simulation, QueryCriteria



# SVET - this is not needed
# we have add_report_simulation_stats in builtin reporters, but it needs a task and manifest
# I'm guessing wherever it is used, it'll have to be moved
class ReportSimulationStats(BaseReport):
    def __init__(self,
                 type="ReportSimulationStats"):
        BaseReport.__init__(self, type=type)

    def to_dict(self):
        d = super(ReportSimulationStats, self).to_dict()
        # d.update({"Report_File_Name": 'ReportSimulationStats' + '.json'})
        return d


def update_basic_params(config, manifest, project_path):
    # SVET - logLevel params might have issues
    vector_species = ['arabiensis', 'funestus', 'gambiae']
    add_species(config, manifest, vector_species)
    config.parameters.Birth_Rate_Dependence = "FIXED_BIRTH_RATE"
    config.parameters.Age_Initialization_Distribution_Type = 'DISTRIBUTION_COMPLEX'
    # SVET - this is how we add parameters that aren't in the schema for sft tests, hope this works
    config.parameters['logLevel_JsonConfigurable'] = 'ERROR'
    config.parameters['logLevel_VectorHabitat'] = 'ERROR'
    config.parameters['logLevel_StandardEventCoordinator'] = 'ERROR'
    config.parameters.Enable_Default_Reporting = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Property_Output = 0

    config.parameters.Enable_Demographics_Risk = 1
    config.parameters.Report_Detection_Threshold_Blood_Smear_Parasites = 0,
    config.parameters.Parasite_Smear_Sensitivity = 0.01  # number of microliters of blood examined
    config.parameters.RDT_Sensitivity = 0.1

    config.parameters.Incubation_Period_Distribution = 'CONSTANT_DISTRIBUTION',
    config.parameters.Incubation_Period_Constant = 3  # parameter shortened from 7 to yield a 12-13 day incubation period
    config.parameters.Immunity_Acquisition_Factor = 1
    config.parameters.Immunity_Initialization_Distribution_Type = "DISTRIBUTION_OFF"
    config.parameters.Immunity_Mortality_Factor = 1
    config.parameters.Immunity_Transmission_Factor = 1

    # read in and set vector bionomics
    vector_bionomics = pd.read_csv(os.path.join(project_path, 'simulation_inputs', 'vector_bionomics.csv'))
    for vector in vector_species:
        bionomics_row = vector_bionomics[vector_bionomics['species'] == vector].reset_index()
        # overwrite=True is only needed for parameters that are lists (as you can either overwrite or add to them)
        set_species_param(config, vector, 'Anthropophily', bionomics_row['Anthropophily'][0])
        set_species_param(config, vector, 'Indoor_Feeding_Fraction', bionomics_row['Indoor_Feeding_Fraction'][0])
        set_species_param(config, vector, 'Days_Between_Feeds', bionomics_row['Days_Between_Feeds'][0])


def habitat_scales(project_path):
    # done
    rel_abundance_fname = os.path.join(project_path, 'ento', 'DS_vector_rel_abundance.csv')
    rdf = pd.read_csv(rel_abundance_fname)
    rdf = rdf.rename(columns={'Anopheles_arabiensis': 'arabiensis_scale_factor',
                              'Anopheles_coluzzii_gambiae': 'gambiae_scale_factor',
                              'Anopheles_funestus_subgroup': 'funestus_scale_factor'})
    rdf = rdf.set_index('DS')
    return rdf


def set_input_files(config, hfca, archetype_hfca, population_size=1000):
    # SVET - might not be able to add these to config
    # i think this is how to add to not purge : config.parameters['Python_Script_Path'] = 'LOCAL'
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters["DS"] = hfca
    config.parameters["Archetype"] = archetype_hfca
    # SVET - demographics filenames will be loaded elsewhere. in the demographics section?

    # SVET - keeping for reference
    # cb.update_params( {
    #     'DS': hfca,
    #     'Archeype': archetype_hfca,
    #     # 'Demographics_Filenames': [os.path.join('demographics_and_climate', hfca, '%s_demographics_wVaxSMC_IPTi.json' % hfca)],
    #     # 'Demographics_Filenames': [os.path.join('demographics_and_climate', hfca, '%s_demographics.json' % hfca)],
    #     'Demographics_Filenames': [os.path.join('demographics_and_climate', '_entire_country', 'demographics_each_admin_%i.json' % population_size)],  # 'season_calib_deomgraphics.json')],
    #     "Climate_Model": "CLIMATE_CONSTANT",
    #     # "Air_Temperature_Filename": os.path.join('demographics_and_climate',archetype_hfca, '%s_air_temperature_daily_2016.bin' % archetype_hfca),
    #     # "Land_Temperature_Filename": os.path.join('demographics_and_climate',archetype_hfca, '%s_air_temperature_daily_2016.bin' % archetype_hfca),
    #     # "Rainfall_Filename": os.path.join('demographics_and_climate',archetype_hfca, '%s_rainfall_daily_2016.bin' % archetype_hfca),
    #     # "Relative_Humidity_Filename": os.path.join('demographics_and_climate',archetype_hfca, '%s_relative_humidity_daily_2016.bin' % archetype_hfca)
    # })

    return {'DS': hfca}


def get_burnin_exp(burnin_id='', burnin_fname=''):
    # SVET - Zhaowei help!
    COMPS_login('https://comps.idmod.org')
    if burnin_id:
        print("building from pickup: identifying burnin addresses from %s" % burnin_id)
        # expt = retrieve_experiment(burnin_id, force_update=True)  # use this if the experiment address is not found even though the experiment succeeded on COMPS
        expt = retrieve_experiment(burnin_id)
        # creating data with all the simulation tags
        ser_df = pd.DataFrame([x.tags for x in expt.simulations])
        # getting paths for all the sims
        # ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])  # old version that no longer works when initial simulation attempts fail
        ser_df["outpath"] = pd.Series(
            [Simulation.get(sim.id, QueryCriteria().select_children('hpc_jobs')).hpc_jobs[-1].working_directory for
             sim in expt.simulations])
    else:
        print("building from pickup: identifying burnin addresses from %s" % burnin_fname)
        ser_df = pd.read_csv(burnin_fname)
    return ser_df


def set_up_hfca(config, manifest, campaign, hfca, archetype_hfca=None,
                pull_from_serialization=False, ser_date=50 * 365,
                hdf=None, lhdf=None, population_size=1000,
                hab_multiplier=-1, run_number=-1, use_arch_burnin=False, ser_df=pd.DataFrame()):
    # SVET - function signature changed
    # this function needs a campaign (for the intervention), config (to set parameters),
    # and manifest(to configure the linear spline in set_habitats())
    # to do all the things it's trying to do.
    # It will probably have to be split up since we separate campaign and config stuff?
    set_input_files(config, hfca, archetype_hfca, population_size)
    if not archetype_hfca:
        archetype_hfca = hfca

    set_habitats(config, manifest, hfca, hdf, lhdf, archetype_hfca, abs(hab_multiplier))

    # SVET - ask Monique
    # I'm not sure what this is trying to do?? There's a "Larval_Habitat_Multiplier" in a
    # "ScaleLarvalHabitat" intervention. I'm not sure what this is doing in config?
    # is this x_Temporary_Larval_Habitat???
    # From docs: LarvalHabitatMultiplier This is read but not used. Both the LarvalHabitatMultiplier and the
    # x_Temporary_Larval_Habitat are applied when a habitat is created at the beginning of a simulation. When the
    # habitat is serialized, it is stored with the results of these multipliers. If you are reading from a serialized
    # file and Serialization_Mask_Node_Read = 0, then you ignore both LarvalHabitatMultiplier and
    # x_Temporary_Larval_Habitat settings and just use what was stored in the serialized file.
    # If Serialization_Mask_Node_Read = 16, then we create new habitats and ignore what is in the serialized file.
    # We do use the x_Temporary_Larval_Habitat setting to adjust the habitat, but it is a known issue that we don’t
    # also use LarvalHabitatMultiplier.
    if pull_from_serialization:
        hab_scale_factor_param_name = 'Habitat_Multiplier'
        if use_arch_burnin:
            if hab_multiplier >= 0 and run_number >= 0:
                ser_df[hab_scale_factor_param_name] = ser_df[hab_scale_factor_param_name].apply(
                    lambda x: np.round(x, 5))
                sdf = ser_df[(ser_df[hab_scale_factor_param_name] >= (np.round(hab_multiplier, 5) - 0.00001)) &
                             (ser_df[hab_scale_factor_param_name] <= (np.round(hab_multiplier, 5) + 0.00001)) &
                             (ser_df['Run_Number'] == run_number) &
                             (ser_df['admin_name'] == archetype_hfca)]
                cb.update_params({
                    hab_scale_factor_param_name: hab_multiplier
                })
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

    # SVET - this will probably need to be pulled out
    add_outbreak_individual(campaign, demographic_coverage=0.002, start_day=35, repetitions=-1,
                            timesteps_between_repetitions=73)

    return {'admin_name': hfca}


def set_habitats(config, manifest, hfca, hdf, lhdf, archetype_hfca, hab_multiplier):
    # done
    a = hdf.at[hfca, 'arabiensis_scale_factor']
    f = hdf.at[hfca, 'funestus_scale_factor']
    g = hdf.at[hfca, 'gambiae_scale_factor']
    tot = a + f + g
    a /= tot
    f /= tot
    g /= tot
    fraction = (max(0.00001, a), max(0.00001, f), max(0.00001, g))

    my_spline, maxvalue, const = load_spline_and_scale_factors(lhdf, archetype_hfca)
    const_mult = 1 if hab_multiplier >= 1 else hab_multiplier

    for (s, sp) in zip(fraction, ['arabiensis', 'funestus', 'gambiae']):
        linear_spline_habitat = configure_linear_spline(manifest,
                                                        max_larval_capacity=pow(10, maxvalue) * s * hab_multiplier,
                                                        capacity_distribution_number_of_years=1,
                                                        capacity_distribution_over_time={
                                                            "Times": [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
                                                                      334],
                                                            "Values": my_spline
                                                        }
                                                        )
        set_species_param(config, sp, "Habitats", linear_spline_habitat)
        # this function updates EMOD parameters to what is requested based on the calibration parameter sampling
        set_species_param(config, sp, "Habitats", {"Habitat_Type": "CONSTANT",
                                                   "Max_Larval_Capacity": pow(10, const) * s * const_mult})


def load_spline_and_scale_factors(lhdf, archetype_hfca):
    # done
    lhdf = lhdf.set_index('archetype')
    my_spline = [lhdf.at[archetype_hfca, 'MonthVal%d' % x] for x in range(1, 13)]
    maxvalue = lhdf.at[archetype_hfca, 'MaxHab']
    const = lhdf.at[archetype_hfca, 'Constant']
    # pop_scale = lhdf.at[archetype_hfca, 'pop_scale']

    return my_spline, maxvalue, const  # , pop_scale


def load_master_csv(project_path):
    # done
    master_csv = os.path.join(project_path, 'admin_pop_archetype.csv')
    df = pd.read_csv(master_csv, encoding='latin')
    df['admin_name'] = df['admin_name'].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    df = df.set_index('admin_name')
    return df


# def load_pop_scale_factor(lhdf, archetype_hfca) :
#
#     df = load_master_csv()
#     my_spline, maxvalue, const, pop_scale = load_spline_and_scale_factors(lhdf, archetype_hfca)
#     scale_factor = 1 / 10000. * (1 / pop_scale)
#     pop = df.at[archetype_hfca, 'population']
#
#     return pop, scale_factor


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


if __name__ == "__main__":
    pass
