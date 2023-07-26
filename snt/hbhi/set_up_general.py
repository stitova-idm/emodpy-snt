import copy
from logging import exception
import os
import pandas as pd
import numpy as np
from dtk.interventions.outbreakindividual import recurring_outbreak
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import update_species_param, set_species, set_larval_habitat
from malaria.interventions.malaria_drugs import set_drug_param
from malaria.reports.MalariaReport import add_filtered_report
from simtools.Utilities.Experiments import retrieve_experiment

def initialize_cb (years, serialize, defaults='MALARIA_SIM', yr_plusone=True,
                   ser_time_step = None, event_reporter = False,
                   filtered_report = None, x_pop_scale=1):
    """Initialize config builder with preset params

    Start a (default) MALARIA_SIM config builder with presets including:
    durations, logging, demographics, serialization, vector and reports.

    Parameters
    ----------
    years: int
        Simulation duration in years
    serialize: bool
        To serialize simulation outcome or not
    defaults: str, default: 'MALARIA_SIM'
        No idea if it works for other defaults...
    yr_plusone: bool, default: True
        For serialization, simulation duration needs to plus one day to work
    ser_time_step: list of int, default: None
        A list defining the sim day(s) to serialize the population. When set
        to None, serialization will occur at `years*365`
    event_reporter: bool, default = False
        If True, switch on 'Report_Event_Recorder' with default event of
        'Received_Severe_Treatment'
    filtered_report: int, default = None
        Number of years (counting from last) to add filtered report on. None
        means no filtered report added
    x_pop_scale: int, default = 1
        Scaling factor for simulation population,
        applies to x_Base_Population, x_Birth, x_Temporary_Larval_Habitat

    Returns
    -------
    DTKConfigBuilder
    
    """
    cb = DTKConfigBuilder.from_defaults(defaults)

    # Run time
    cb.update_params({'Simulation_Duration': years*365+yr_plusone})
    
    # Logging
    cb.update_params({
        'logLevel_JsonConfigurable' : 'ERROR',
        'logLevel_VectorHabitat' : 'ERROR',
        'logLevel_StandardEventCoordinator' : 'ERROR',
        'logLevel_SusceptibilityMalaria' : 'ERROR'
    })

    # Demographics
    cb.update_params({
        "Birth_Rate_Dependence": "FIXED_BIRTH_RATE",
        "Age_Initialization_Distribution_Type": 'DISTRIBUTION_COMPLEX',
        'Disable_IP_Whitelist' : 1,
        'x_Base_Population': x_pop_scale,
        'x_Birth': x_pop_scale
    }) 

    # Serialization
    if ser_time_step is None:
       ser_time_step = [365*years] 
    if serialize :
        cb.update_params({
            'Serialization_Time_Steps' : ser_time_step,
            'Serialization_Type': 'TIMESTEP',
            'Serialized_Population_Writing_Type' : 'TIMESTEP',
            'Serialization_Mask_Node_Write': 0,  
            # 0 corresponds to the previous version default: the same larval habitat 
            # parameters will be used in the burnin and pickup (from the burnin config)
            'Serialization_Precision': 'REDUCED'
        })
    else:
        cb.update_params({
            'Serialization_Type': 'NONE',
            'Serialized_Population_Writing_Type': 'NONE'
        })
    
    # Vector
    cb.update_params({
        "Vector_Species_Names": ['arabiensis', 'funestus', 'gambiae'],
        'x_Temporary_Larval_Habitat': 0.2 * x_pop_scale # Set as 0.2 for GN and NG compatibility, future work should be 1
    })
    set_species(cb, ['arabiensis', 'funestus', 'gambiae'])

    update_species_param(cb, 'arabiensis', 'Anthropophily', 0.65, overwrite=True)
    update_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
    update_species_param(cb, 'funestus', 'Anthropophily', 0.35, overwrite=True)
    update_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 0.92, overwrite=True)
    update_species_param(cb, 'gambiae', 'Anthropophily', 0.85, overwrite=True)
    update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.95, overwrite=True)

    # Report
    cb.update_params( {
        'Enable_Default_Reporting': 0,
        'Enable_Demographics_Risk': 1,
        'Enable_Property_Output': 0,
        'Enable_Vector_Species_Report': 0,
        'Report_Detection_Threshold_Blood_Smear_Parasites': 50,
        "Parasite_Smear_Sensitivity": 0.02,  # 50/uL
        'RDT_Sensitivity': 0.1
    })

    if event_reporter:
        cb.update_params({
            "Report_Event_Recorder": 1,
            "Report_Event_Recorder_Individual_Properties": [],
            "Report_Event_Recorder_Events" : ['Received_Severe_Treatment',
                                              # 'Received_Treatment', 'NewClinicalCase', 
                                              # 'NewSevereCase', 'Received_Campaign_Drugs', 
                                              # 'No_SMC_Fever'
                                              ],
            "Report_Event_Recorder_Ignore_Events_In_List": 0
        })
    
    if filtered_report:
        num_year = filtered_report
        start = yr_plusone + (years - num_year) * 365
        end = yr_plusone + years * 365
        add_filtered_report(cb, start=start, end=end)
    
    return cb

def load_master_csv(projectpath, file = None, country = None) :
    if file is None:
        if country in ['burkina', 'guinea']:
            fname = country + "_DS_pop.csv"
            DS_Name = 'DS_Name'
        elif country == 'nigeria':
            fname = country + "_LGA_pop.csv"
            DS_Name = 'LGA'
        elif country is None:
            raise Exception("Must specify either the file or the country.")
        else:
           raise Exception("Only burkina, guinea and nigeria are supported for country specification.")
    else:
        fname = file
        DS_Name = 'DS_Name'

    master_csv = os.path.join(projectpath, fname)
    df = pd.read_csv(master_csv, encoding='latin')
    df['DS_Name'] = df[DS_Name].str.normalize('NFKD').str.encode('ascii', errors='ignore').str.decode('utf-8')
    df = df.set_index('DS_Name')
    return df

def set_input_files(cb, my_ds, archetype_ds=None, demographic_suffix='', 
                    climate_suffix='', climate_prefix=True, use_archetype=True) :
    if archetype_ds is None:
        archetype_ds = my_ds

    if demographic_suffix is not None:
        if not demographic_suffix.startswith('_') and not demographic_suffix == '':
            demographic_suffix = '_' + demographic_suffix
    
    if climate_suffix is not None:
        if not climate_suffix.startswith('_') and not climate_suffix == '':
            climate_suffix = '_' + climate_suffix
   
    if use_archetype:
        ds = archetype_ds
    else:
        ds = my_ds

    cb.update_params({
        'District_Sanitaire' : my_ds,
        'Archetype' : archetype_ds
    })
    if demographic_suffix is not None:
        cb.update_params({
            'Demographics_Filenames': [os.path.join(ds, '%s_demographics%s.json' % (ds, demographic_suffix))]
        })
    if climate_suffix is not None:
        if climate_prefix:
            cb.update_params({
                "Air_Temperature_Filename": os.path.join(ds, '%s_air_temperature_daily%s.bin' % (ds, climate_suffix)),
                "Land_Temperature_Filename": os.path.join(ds, '%s_air_temperature_daily%s.bin' % (ds, climate_suffix)),
                "Rainfall_Filename": os.path.join(ds, '%s_rainfall_daily%s.bin' % (ds, climate_suffix)),
                "Relative_Humidity_Filename": os.path.join(ds, '%s_relative_humidity_daily%s.bin' % (ds, climate_suffix))
            })
        else:
            cb.update_params({
                "Air_Temperature_Filename": os.path.join(ds, 'air_temperature_daily%s.bin' % climate_suffix),
                "Land_Temperature_Filename": os.path.join(ds, 'air_temperature_daily%s.bin' % climate_suffix),
                "Rainfall_Filename": os.path.join(ds, 'rainfall_daily%s.bin' % climate_suffix),
                "Relative_Humidity_Filename": os.path.join(ds, 'relative_humidity_daily%s.bin' % climate_suffix)
            })

    return {'DS_Name' : my_ds}

def setup_ds(cb, my_ds, archetype_ds=None,
             pull_from_serialization=False,
             burnin_id='', ser_date=50*365,
             burnin_fname='',
             rel_abund_df=None, lhdf=None, use_arch_burnin=True,
             from_arch=None, demographic_suffix='',
             climate_suffix='',
             climate_prefix=True,
             use_arch_input=True,
             hab_multiplier=1, run_number=-1,
             parser_default='HPC',DS_Name='DS_Name',
             serialize_match_tag=None,
             serialize_match_val=None) :
    """Setting an individual DS up to add to config builder

    Given a DS, pull its archetype, demographic, climate and vector files, 
    its corresponding serialized population (if applicable) and other info
    to feed into the config builder. Usually used within `ModFn` and 
    `ModBuilder`.

    Parameters
    ----------
    cb : DTKConfigBuilder
    my_ds : str
    archetype_ds : str, default: None
        The archetype DS of `my_ds`
    pull_from_serialization : bool, default: False
    burnin_id : str, default: ''
        Experiment ID to pull serialization from
    ser_date : int, default: 50*365
        Number of sim days in the serialized population
    burnin_fname : str, default: ''
        Path to csv files containing dataframe that specifies where to look
        up for serialized populations; Only needed if `pull_from_serialization`
        is True and ignore if `burnin_id` is specified
    rel_abund_df : pandas.DataFrame, default: None
        Data frame of relative vector abundance by DS (i.e., interspecies 
        relative abundance)
    lhdf : pandas.DataFrame, default: None
        Data frame of larval habitats per month per DS
    use_arch_burnin : bool, default: True
        Pull serialization from corresponding archetype instead of DS itself
    from_arch : bool, default: None
        (deprecated) Same as `use_arch_burnin`, supersedes its value if not
        `None`
    demographic_suffix : str, default: ''
        Default demographic input file is '(DSNAME)_demographic', this argument
        adds suffix to input file name. Note: if set to `None`, demographic input
        file will not be added
    climate_suffix : str, default: ''
        Same behaviour as `demographic_suffix`
    climate_prefix : bool, default: True
        Whether to add DS name in front of climate files.
    use_arch_input : bool, default: True
        If True, use archetype DS's demographic and climate input file instead
    hab_multiplier : float, default: 1
    run_number : int, default: -1
        `run_number` will be matched with `run_number` of the serialized population
        from the burnin experiment
    parser_default : {'HPC', 'NUCLUSTER', 'LOCAL'}
        Specifies where does the serialized population live
    DS_Name : str, default: 'DS_Name'
        The variable name that tags a simulation run to a DS; Could be 'LGA' in the
        case of Nigeria
    
    Returns
    -------
    dict
        Dictionary of two keys, the DS_Name and the corresponding archetype DS
    
    """

    # For backward compatibility, from_arch supersedes use_arch_burnin
    if from_arch is not None:
        use_arch_burnin = from_arch
    set_input_files(cb, my_ds, archetype_ds, demographic_suffix, climate_suffix, climate_prefix, use_arch_input)
    if not archetype_ds :
        archetype_ds = my_ds
    
    if rel_abund_df is not None:
        set_habitats(cb, rel_abund_df, lhdf, archetype_ds, abs(hab_multiplier)) # To be fixed

    if pull_from_serialization:
        # serialize_match check
        if serialize_match_tag:
            if serialize_match_val:
                if not len(serialize_match_tag) == len(serialize_match_val):
                    raise Exception("serialize_match_tag and serialize_match_val not same length.")
            else:
                raise Exception("serialize_match_val must be specified if serialize_match_tag is not None.")
        else:
            serialize_match_tag = ['Habitat_Multiplier', 'Run_Number']
            serialize_match_val = [hab_multiplier, run_number]

        if parser_default == 'HPC':
            from simtools.Utilities.COMPSUtilities import COMPS_login
            COMPS_login('https://comps.idmod.org')
        # print("building from pickup")
        
        # serialization
        # print("retrieving burnin")
        if burnin_id :
            expt = retrieve_experiment(burnin_id)
            # creating data with all the simulation tags
            ser_df = pd.DataFrame([x.tags for x in expt.simulations])
            # getting paths for all the sims
            ser_df["outpath"] = pd.Series([sim.get_path() for sim in expt.simulations])
        else :
            ser_df = pd.read_csv(burnin_fname)

        if use_arch_burnin :
            ser_df = ser_df[ser_df[DS_Name] == archetype_ds]
        else :
            ser_df = ser_df[ser_df[DS_Name] == my_ds]
        
        sdf = ser_df.copy()
        for t, v in zip(serialize_match_tag, serialize_match_val):
            if type(v) == float:
                sdf[t] = sdf[t].apply(lambda x : np.round(x, 5))
                sdf = sdf[sdf[t] == np.round(v, 5)]
            else:
                sdf = sdf[sdf[t] == v]
            cb.update_params({t: v})
        ser_path = sdf['outpath'].values[0]

        cb.update_params( {
            'Serialized_Population_Reading_Type': 'READ',
            'Serialized_Population_Path' : os.path.join(ser_path, 'output'),
            'Serialized_Population_Filenames' : ['state-%05d.dtk' % ser_date],
            'Enable_Random_Generator_From_Serialized_Population' : 0,
            'Serialization_Mask_Node_Read' : 0
        })
    else:
        cb.update_params({
            'Serialized_Population_Reading_Type': 'NONE',
        })
    
    recurring_outbreak(cb, start_day=182, outbreak_fraction=0.01, tsteps_btwn=365)

    return {DS_Name : my_ds,
            'archetype' : archetype_ds}


def set_habitats(cb, hdf, lhdf, archetype_ds, hab_multiplier, my_ds=None) :
    if my_ds is not None:
        ds = my_ds
    else:
        ds = archetype_ds
    hdf1 = hdf.loc[ds]
    scale_factors = ['arabiensis_scale_factor', 'funestus_scale_factor', 'gambiae_scale_factor']
    prop_list = list(hdf1[scale_factors]/np.sum(hdf1[scale_factors]))
    fraction = tuple(prop_list)

    ls_hab_ref = {'Capacity_Distribution_Number_Of_Years': 1,
                  'Capacity_Distribution_Over_Time': {
                      'Times': [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334],
                      'Values': []
                  },
                  'Max_Larval_Capacity': 800000000}

    my_spline, maxvalue, const = load_spline_and_scale_factors(lhdf, archetype_ds)
    const_mult = 1 if hab_multiplier >=1 else hab_multiplier

    for (s, sp) in zip(fraction, ['arabiensis', 'funestus', 'gambiae']) :
        hab = copy.copy(ls_hab_ref)
        hab['Capacity_Distribution_Over_Time']['Values'] = my_spline
        hab['Max_Larval_Capacity'] = pow(10, maxvalue)*s*hab_multiplier
        # this function updates EMOD parameters to what is requested based on the calibration parameter sampling
        set_larval_habitat(cb, { sp : {'LINEAR_SPLINE' : hab,
                                       'CONSTANT' : pow(10, const)*s*const_mult}})


def load_spline_and_scale_factors(lhdf, archetype_ds) :

    if lhdf is None :
        my_spline = [0]*12
        maxvalue = 1
        const = 1

    else :
        lhdf = lhdf.set_index('archetype')
        my_spline = [lhdf.at[archetype_ds, 'Month%d' % x] for x in range(1, 13)]
        maxvalue = lhdf.at[archetype_ds, 'MaxHab']
        const = lhdf.at[archetype_ds, 'Constant']

    return my_spline, maxvalue, const

# Deprecated, included here as legacy. 
# Input DS_vector_rel_abundance should standardize to the format of GN
def habitat_scales(mdf, projectpath, write=False) :
    df = mdf.reset_index()
    try:
        rel_abundance_fname = os.path.join(projectpath, 'SpatialClustering_BF', 'DS_vector_rel_abundance.csv')
    except:
        rel_abundance_fname = os.path.join(projectpath, 'SpatialClustering', 'DS_vector_rel_abundance.csv')
    rdf = pd.read_csv(rel_abundance_fname)
    del rdf['DS_Name']
    rdf = pd.merge(rdf, df[['DS', 'DS_Name']], on='DS')
    rdf = rdf.rename(columns={'Anopheles_arabiensis' : 'arabiensis_scale_factor',
                              'Anopheles_coluzzii_gambiae' : 'gambiae_scale_factor',
                              'Anopheles_funestus_subgroup' : 'funestus_scale_factor'})
    rdf = rdf.set_index('DS_Name')
    if write:
        rdf.to_csv(os.path.join(projectpath, 'simulation_inputs', 'DS_vector_rel_abundance.csv'))
    return rdf


def load_rel_abund_df(projectpath) :

    rel_abundance_fname = os.path.join(projectpath, 'simulation_inputs', 'DS_vector_rel_abundance.csv')
    rel_abund_df = pd.read_csv(rel_abundance_fname)
    rel_abund_df = rel_abund_df.set_index('DS_Name')
    return rel_abund_df

def set_spaq_params(cb) :

    set_drug_param(cb, 'Amodiaquine', "Drug_PKPD_C50", 1)
    set_drug_param(cb, 'Amodiaquine', "Max_Drug_IRBC_Kill", 0.23)
    set_drug_param(cb, 'Sulfadoxine', "Drug_PKPD_C50", 0.2*6)
    set_drug_param(cb, 'Sulfadoxine', "Max_Drug_IRBC_Kill", 0.506*0.675)
    set_drug_param(cb, 'Pyrimethamine', "Drug_PKPD_C50", 8*6)
    set_drug_param(cb, 'Pyrimethamine', "Max_Drug_IRBC_Kill", 0.6*0.6417)

    set_drug_param(cb, 'Amodiaquine', "Drug_Cmax", 95)
    set_drug_param(cb, 'Amodiaquine', "Drug_Decay_T1", 0.775)
    set_drug_param(cb, 'Amodiaquine', "Drug_Decay_T2", 37.5)
    set_drug_param(cb, 'Amodiaquine', "Drug_Vd", 15.6)
    set_drug_param(cb, 'Amodiaquine', "Fractional_Dose_By_Upper_Age", [{"Upper_Age_In_Years": 1, "Fraction_Of_Adult_Dose": 0.376}])

    return {'Drug': "User-defined (set_spaq_params)"}
