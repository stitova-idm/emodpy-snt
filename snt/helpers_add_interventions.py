import os
import copy
import warnings
import pandas as pd
import numpy as np
from datetime import date
from malaria.interventions.health_seeking import add_health_seeking
from dtk.interventions.outbreakindividual import recurring_outbreak
from dtk.vector.species import update_species_param
from dtk.vector.species import set_larval_habitat
from simtools.Utilities.Experiments import retrieve_experiment
from simtools.Utilities.COMPSUtilities import COMPS_login
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.biting_risk import change_biting_risk
from dtk.interventions.irs import add_IRS
from dtk.interventions.property_change import change_individual_property

# from dtk.interventions.habitat_scale import scale_larval_habitats
from malaria.reports.MalariaReport import add_event_counter_report
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign, add_diagnostic_survey
from malaria.interventions.malaria_vaccdrug_campaigns import add_vaccdrug_campaign
from malaria.interventions.adherent_drug import configure_adherent_drug
from dtk.interventions.property_change import change_individual_property_at_age
from simulation.helpers_sim_setup import update_smc_access_ips
from simulation.load_paths import load_box_paths
from malaria.interventions.malaria_vaccine import add_vaccine


def add_hfca_hs(cb, hs_df, hfca, seed_index=0):

    #df = hs_df[hs_df['repDS'] == hfca]
    if 'LGA' in hs_df.columns.values:
        df = hs_df[hs_df['LGA'] == hfca]
    elif 'repDS' in hs_df.columns.values:
        df = hs_df[hs_df['repDS'] == hfca]
    else:
        df = hs_df[hs_df['admin_name'] == hfca]
    if 'seed' in df.columns.values:
        df = df[df['seed'] == seed_index]
    for r, row in df.iterrows() :
        add_hs_from_file(cb, row)

    return len(df)


def add_hs_from_file(cb, row):

    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    start_day = row['simday']
    severe_cases = row['severe_coverage']

    add_health_seeking(cb, start_day=start_day,
                       targets=[{'trigger': 'NewClinicalCase', 'coverage': hs_child, 'agemin': 0, 'agemax': 5,
                                 'seek': 1, 'rate': 0.3},
                                {'trigger': 'NewClinicalCase', 'coverage': hs_adult, 'agemin': 5, 'agemax': 100,
                                 'seek': 1, 'rate': 0.3},
                                ],
                       drug=['Artemether', 'Lumefantrine'], duration=row['duration'])
    add_health_seeking(cb, start_day=start_day,
                       targets=[{'trigger': 'NewSevereCase', 'coverage': severe_cases, 'seek': 1, 'rate': 0.5}], #change by adding column and reviewing literature
                       drug=['Artemether', 'Lumefantrine'], duration=row['duration'],
                       broadcast_event_name = 'Received_Severe_Treatment')

def add_nmf_hs(cb, hs_df, nmf_df, hfca, seed_index=0):

    # if no NMF rate is specified, assume all age groups have 0.0038 probability each day
    if nmf_df.empty:
        nmf_df = pd.DataFrame({'U5_nmf': [0.0038], 'adult_nmf': [0.0038]})
    elif nmf_df.shape[0] != 1:
        warnings.warn('The NMF dataframe has more than one row. Only values in the first row will be used.')
    nmf_row = nmf_df.iloc[0]


    # apply the health-seeking rate for clinical malaria to NMFs
    df = hs_df[hs_df['admin_name'] == hfca]
    if 'seed' in df.columns.values:
        df = df[df['seed'] == seed_index]
    for r, row in df.iterrows():
        add_nmf_hs_from_file(cb, row, nmf_row)

    return len(df)


def add_nmf_hs_from_file(cb, row, nmf_row):

    hs_child = row['U5_coverage']
    hs_adult = row['adult_coverage']
    start_day = row['simday']
    duration = row['duration']
    if start_day == 0:  # due to dtk diagnosis/treatment configuration, a start day of 0 is not supported
        start_day = 1  # start looking for NMFs on day 1 (not day 0) of simulation
        if duration > 1:
            duration = duration - 1
    nmf_child = nmf_row['U5_nmf']
    nmf_adult = nmf_row['adult_nmf']
    if duration < 1000:
        start_day_list = [start_day]
        duration_list = [duration]
    else:
        interval = 900  # should be less than the EMOD max (which is currently 1000)
        start_day_list_0 = list(range(0, duration, interval))
        start_day_list = [start_day + ss for ss in start_day_list_0]
        duration_list = [interval if ((ss + interval) < duration) else (duration - ss) for ss in start_day_list_0]

    for ii in range(len(start_day_list)):
        add_drug_campaign(cb, 'MSAT', 'AL', start_days=[start_day_list[ii]],
                          target_group={'agemin': 0, 'agemax': 5},
                          coverage=nmf_child * hs_child,
                          repetitions=duration_list[ii], tsteps_btwn_repetitions=1,
                          diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                          receiving_drugs_event_name='Received_NMF_Treatment')
        add_drug_campaign(cb, 'MSAT', 'AL', start_days=[start_day_list[ii]],
                          target_group={'agemin': 5, 'agemax': 120},
                          coverage=nmf_adult * hs_adult,
                          repetitions=duration_list[ii], tsteps_btwn_repetitions=1,
                          diagnostic_type='PF_HRP2', diagnostic_threshold=5,
                          receiving_drugs_event_name='Received_NMF_Treatment')


def add_hfca_irs(cb, irs_df, hfca, seed_index=0) :
    irs_df = irs_df[irs_df['admin_name'].str.upper() == hfca.upper()]
    if 'seed' in irs_df.columns.values:
        irs_df = irs_df[irs_df['seed'] == seed_index]
    for r, row in irs_df.iterrows() :
        add_IRS(cb, start=row['simday'], coverage_by_ages=[{"coverage": row['effective_coverage'], "min": 0, "max": 100}],
                killing_config={
                    "class": "WaningEffectExponential",
                    "Decay_Time_Constant": row['mean_duration'],
                    "Initial_Effect": row['initial_kill']},
                )

    return len(irs_df)

# all itn
def add_hfca_itns(cb, itn_df, itn_anc_df, itn_anc_adult_birthday_years, itn_epi_df, itn_chw_df, itn_chw_annual_df, hfca, itn_use_seasonality, itn_decay_params, seed_index=0) :#
    if not itn_df.empty:
        df = itn_df[itn_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets = len(df)
        for r, row in df.iterrows() :
            add_itn_from_file(cb, row, itn_use_seasonality, itn_decay_params)
    else:
        nets = 0

    if not itn_anc_df.empty:
        df = itn_anc_df[itn_anc_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        add_itn_anc(cb, df, itn_anc_adult_birthday_years, itn_use_seasonality, itn_decay_params)
        nets += len(df)

    if not itn_epi_df.empty:
        df = itn_epi_df[itn_epi_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        add_birthday_routine_itn_from_file(cb, df, itn_use_seasonality, itn_decay_params)
        nets += len(df)

    if not itn_chw_df.empty:
        df = itn_chw_df[itn_chw_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets += len(df)
        for r, row in df.iterrows() :
            add_monthly_CHW_dist(cb, row, itn_use_seasonality, itn_decay_params)

    if not itn_chw_annual_df.empty:
        df = itn_chw_annual_df[itn_chw_annual_df['admin_name'].str.upper() == hfca.upper()]
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        df = df.drop_duplicates()
        nets += len(df)
        for r, row in df.iterrows() :
            add_annual_CHW_dist(cb, row, itn_use_seasonality, itn_decay_params)

    return nets


def add_itn_from_file(cb, row, itn_use_seasonality, itn_decay_params):
    # net retention time, decay rates of killing and blocking
    itn_Expiration_Period_Log_Normal_Mu = row['net_life_lognormal_mu']
    itn_Expiration_Period_Log_Normal_Sigma = row['net_life_lognormal_sigma']
    itn_decay_kill = itn_decay_params['kill_decay_time'][0]
    itn_decay_block = itn_decay_params['block_decay_time'][0]

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']
    seasonal_offset = row['simday'] % 365
    seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

    zipped_lists = zip(seasonal_times, seasonal_scales)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
    if seasonal_times[0] > 0:
        seasonal_times.insert(0, 0)
        seasonal_scales.insert(0, seasonal_scales[-1])

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    if coverage_all > 0:
        age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]
    else:
        age_dep = [1 for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    add_ITN_age_season(cb, start=row['simday'], demographic_coverage=coverage_all,
                       killing_config={
                           "Initial_Effect": row['kill_initial'],
                           "Decay_Time_Constant": itn_decay_kill,
                           "class": "WaningEffectExponential"},
                       blocking_config={
                           "Initial_Effect": row['block_initial'],
                           "Decay_Time_Constant": itn_decay_block,
                           "class": "WaningEffectExponential"},
                       discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                              "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                       age_dependence={'Times': [0, 5, 10, 15, 20],
                                       'Values': [x * indoor_net_protection for x in age_dep]},
                       seasonal_dependence={"Times": seasonal_times,
                                            "Values": seasonal_scales}
                       )

    if 'ITN coverage severe treatment' in row :
        if row['simday'] == 895 :
            add_ITN_age_season(cb, start=0, demographic_coverage=row['ITN coverage severe treatment'],
                               killing_config={
                                   "Initial_Effect": row['kill_initial'],
                                   "Decay_Time_Constant": itn_decay_kill,
                                   "class": "WaningEffectExponential"},
                               blocking_config={
                                   "Initial_Effect": row['block_initial'],
                                   "Decay_Time_Constant": itn_decay_block,
                                   "class": "WaningEffectExponential"},
                               discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                              "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                               trigger_condition_list=['Received_Severe_Treatment'],
                               duration=-1,
                               age_dependence={'Times': [0, 5, 10, 15, 20],
                                               'Values': [x * indoor_net_protection for x in age_dep]},
                               seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales}
                               )
    if 'ITN coverage uncomplicated treatment' in row :
        if row['simday'] == 895 :
            add_ITN_age_season(cb, start=0, demographic_coverage=row['ITN coverage uncomplicated treatment'],
                               killing_config={
                                   "Initial_Effect": row['kill_initial'],
                                   "Decay_Time_Constant": itn_decay_kill,
                                   "class": "WaningEffectExponential"},
                               blocking_config={
                                   "Initial_Effect": row['block_initial'],
                                   "Decay_Time_Constant": itn_decay_block,
                                   "class": "WaningEffectExponential"},
                               discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                              "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                               trigger_condition_list=['Received_Treatment'],
                               duration=-1,
                               age_dependence={'Times': [0, 5, 10, 15, 20],
                                               'Values': [x * indoor_net_protection for x in age_dep]},
                               seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales}
                               )


def add_birthday_routine_itn_from_file(cb, itn_epi_df, itn_use_seasonality, itn_decay_params):
    # decay rates of killing and blocking
    itn_decay_kill = itn_decay_params['kill_decay_time'][0]
    itn_decay_block = itn_decay_params['block_decay_time'][0]

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']

    for r, row in itn_epi_df.iterrows() :
        # net retention time
        itn_Expiration_Period_Log_Normal_Mu = row['net_life_lognormal_mu']
        itn_Expiration_Period_Log_Normal_Sigma = row['net_life_lognormal_sigma']
        indoor_net_protection = row['indoor_net_protection']

        # age of individuals receiving nets
        birthday_age = row['birthday_age']

        seasonal_offset = row['simday'] % 365
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        add_ITN_age_season(cb, start=row['simday'], demographic_coverage=row['coverage'],
                           age_min=(birthday_age-0.5), age_max=(birthday_age+0.5),
                           killing_config={
                               "Initial_Effect": row['kill_initial'],
                               "Decay_Time_Constant": itn_decay_kill,
                               "class": "WaningEffectExponential"},
                           blocking_config={
                               "Initial_Effect": row['block_initial'],
                               "Decay_Time_Constant": itn_decay_block,
                               "class": "WaningEffectExponential"},
                           discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                          "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                          "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                           trigger_condition_list=['HappyBirthday'],
                           duration=row['duration'],
                           age_dependence={'Times': [0, 120],
                                           'Values': [indoor_net_protection, indoor_net_protection]},
                           seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales}
                           )


def add_itn_anc(cb, itn_anc_df, itn_anc_adult_birthday_years, itn_use_seasonality, itn_decay_params):

    # decay rates of killing and blocking
    itn_decay_kill = itn_decay_params['kill_decay_time'][0]
    itn_decay_block = itn_decay_params['block_decay_time'][0]

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']

    for r, row in itn_anc_df.iterrows() :
        # net retention time
        itn_Expiration_Period_Log_Normal_Mu = row['net_life_lognormal_mu']
        itn_Expiration_Period_Log_Normal_Sigma = row['net_life_lognormal_sigma']
        indoor_net_protection = row['indoor_net_protection']

        seasonal_offset = row['simday'] % 365
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        # ITNs protecting infants
        add_ITN_age_season(cb, start=row['simday'], demographic_coverage=row['coverage'],
                           killing_config={
                               "Initial_Effect": row['kill_initial'],
                               "Decay_Time_Constant": itn_decay_kill,
                               "class": "WaningEffectExponential"},
                           blocking_config={
                               "Initial_Effect": row['block_initial'],
                               "Decay_Time_Constant": itn_decay_block,
                               "class": "WaningEffectExponential"},
                           discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                          "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                          "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                           birth_triggered=True,
                           duration=row['duration'],
                           age_dependence={'Times': [0, 5],
                                           'Values': [indoor_net_protection, indoor_net_protection]},
                           seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales}
                           )

        # ITNs protecting mothers (since we don't simulate pregnancy explicitly, we approximate ANC ITN coverage among
        #    women by giving them ITNs at certain birthdays, using the ANC ITN coverage for that admin)
        for birthday_year in itn_anc_adult_birthday_years:
            add_ITN_age_season(cb, start=row['simday'], demographic_coverage=row['coverage'],
                               age_min=(birthday_year-0.5), age_max=(birthday_year+0.5), target_gender='Female',
                               killing_config={
                                   "Initial_Effect": row['kill_initial'],
                                   "Decay_Time_Constant": itn_decay_kill,
                                   "class": "WaningEffectExponential"},
                               blocking_config={
                                   "Initial_Effect": row['block_initial'],
                                   "Decay_Time_Constant": itn_decay_block,
                                   "class": "WaningEffectExponential"},
                               discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                              "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                               trigger_condition_list=['HappyBirthday'],
                               duration=row['duration'],
                               age_dependence={'Times': [0, 120],
                                               'Values': [indoor_net_protection, indoor_net_protection]},
                               seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales}
                               )



def add_monthly_CHW_dist(cb, row, itn_use_seasonality, itn_decay_params):
    # net retention time, decay rates of killing and blocking
    itn_Expiration_Period_Log_Normal_Mu = row['net_life_lognormal_mu']
    itn_Expiration_Period_Log_Normal_Sigma = row['net_life_lognormal_sigma']
    itn_decay_kill = itn_decay_params['kill_decay_time'][0]
    itn_decay_block = itn_decay_params['block_decay_time'][0]

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    # add distribution every month for a year (every 30 days starting from simday)
    simday = row['simday']
    for mm in range(12):
        chw_check_day = simday + 30 * mm

        # seasonality in ITN use
        seasonal_scales = itn_use_seasonality['itn_use_scalar']
        seasonal_days = itn_use_seasonality['day']
        seasonal_offset = chw_check_day
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0:
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        add_ITN_age_season(cb, start=chw_check_day, demographic_coverage=coverage_all,
                           killing_config={
                               "Initial_Effect": row['kill_initial'],
                               "Decay_Time_Constant": itn_decay_kill,
                               "class": "WaningEffectExponential"},
                           blocking_config={
                               "Initial_Effect": row['block_initial'],
                               "Decay_Time_Constant": itn_decay_block,
                               "class": "WaningEffectExponential"},
                           discard_times={"Expiration_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                                  "Expiration_Period_Log_Normal_Mu": itn_Expiration_Period_Log_Normal_Mu,
                                                  "Expiration_Period_Log_Normal_Sigma": itn_Expiration_Period_Log_Normal_Sigma},
                           age_dependence={'Times': [0, 5, 10, 15, 20],
                                           'Values': [x * indoor_net_protection for x in age_dep]},
                           seasonal_dependence={"Times": seasonal_times,
                                                "Values": seasonal_scales},
                           dont_allow_duplicates=1
                           )




def add_annual_CHW_dist(cb, row, itn_use_seasonality, itn_decay_params):
    # net retention time, decay rates of killing and blocking
    itn_Expiration_Period_Log_Normal_Mu = row['net_life_lognormal_mu']
    itn_Expiration_Period_Log_Normal_Sigma = row['net_life_lognormal_sigma']
    itn_decay_kill = itn_decay_params['kill_decay_time'][0]
    itn_decay_block = itn_decay_params['block_decay_time'][0]

    # use-coverage by age
    itn_u5 = row['itn_u5']
    itn_5_10 = row['itn_5_10']
    itn_10_15 = row['itn_10_15']
    itn_15_20 = row['itn_15_20']
    itn_o20 = row['itn_o20']
    # a single distribution coverage is given to all age groups, but use differs by age.
    #   set the distribution coverage to be the maximum use coverage across all age groups, then use age_dependence
    #   uses to end up with the appropriate use for that age group
    coverage_all = max([itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20])
    age_dep = [x / coverage_all for x in [itn_u5, itn_5_10, itn_10_15, itn_15_20, itn_o20]]

    # fraction of indoor time protected by net
    indoor_net_protection = row['indoor_net_protection']

    # add distribution every month for a year (every 30 days starting from simday)
    simday = row['simday']

    # seasonality in ITN use
    seasonal_scales = itn_use_seasonality['itn_use_scalar']
    seasonal_days = itn_use_seasonality['day']
    seasonal_offset = simday % 365
    seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in seasonal_days]

    zipped_lists = zip(seasonal_times, seasonal_scales)
    sorted_pairs = sorted(zipped_lists)
    tuples = zip(*sorted_pairs)
    seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
    if seasonal_times[0] > 0:
        seasonal_times.insert(0, 0)
        seasonal_scales.insert(0, seasonal_scales[-1])

    add_ITN_age_season(cb, start=simday, demographic_coverage=coverage_all,
                       killing_config={
                           "Initial_Effect": row['kill_initial'],
                           "Decay_Time_Constant": itn_decay_kill,
                           "class": "WaningEffectExponential"},
                       blocking_config={
                           "Initial_Effect": row['block_initial'],
                           "Decay_Time_Constant": itn_decay_block,
                           "class": "WaningEffectExponential"},
                       discard_times={"Expiration_Period_Distribution": "CONSTANT_DISTRIBUTION",
                                              "Expiration_Period_Constant": 365},
                       age_dependence={'Times': [0, 5, 10, 15, 20],
                                       'Values': [x * indoor_net_protection for x in age_dep]},
                       seasonal_dependence={"Times": seasonal_times,
                                            "Values": seasonal_scales}
                       )


def add_hfca_vaccsmc(cb, smc_df, hfca, effective_coverage_resistance_multiplier=1, seed_index=0):
    df = smc_df[smc_df['admin_name'] == hfca]
    if df.shape[0] > 0:
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]

    if len(df) == 0:
        return len(df)

    for r, row in df.iterrows():
        if 'max_age' in df.columns.values:
            max_smc_age = row['max_age']
        else:
            max_smc_age = 5
        if 'TAT' in df.columns.values:
            TAT = row['TAT']
        else:
            TAT = 0
        if TAT:
            raise ValueError('TAT not yet implemented for vaccSMC')

        # get coverage and ages for SMC in different groups
        cov_high_U5 = row['coverage_high_access_U5'] * effective_coverage_resistance_multiplier
        cov_low_U5 = row['coverage_low_access_U5'] * effective_coverage_resistance_multiplier
        if 'coverage_high_access_5_10' in df.columns.values:
            cov_high_5_10 = row['coverage_high_access_5_10'] * effective_coverage_resistance_multiplier
            cov_low_5_10 = row['coverage_low_access_5_10'] * effective_coverage_resistance_multiplier
            # create list of values for all groups receiving SMC with order: [U5-high, U5-low, O5-high, O5-low]
            age_min_list = [0.25, 0.25, 5, 5]
            age_max_list = [5, 5, max(10, max_smc_age), max(10, max_smc_age)]
            cov_list = [cov_high_U5, cov_low_U5, cov_high_5_10, cov_low_5_10]
            access_list = ['High', 'Low', 'High', 'Low']
        else:
            # create list of values for all groups receiving SMC with order: [U5-high, U5-low, O5-high, O5-low]
            age_min_list = [0.25, 0.25]
            age_max_list = [5, 5]
            cov_list = [cov_high_U5, cov_low_U5]
            access_list = ['High', 'Low']

        for ii in range(len(cov_list)):
            add_vaccdrug_campaign(cb, campaign_type='SMC', start_days=[row['simday']],
                                  coverages=cov_list[ii],
                                  target_group={'agemin': age_min_list[ii],
                                                'agemax': age_max_list[ii]},
                                  ind_property_restrictions=[{'SMCAccess': access_list[ii]}],
                                  receiving_drugs_event=False)  ## If False uses vaccSMC with automatic offset of 17 days, if True, uses vaccDrugSMC

    return len(df)


def smc_adherent_configuration(cb, adherence, sp_resist_day1_multiply):
    smc_adherent_config = configure_adherent_drug(cb,
                                                  doses=[["SulfadoxinePyrimethamine", 'Amodiaquine'],
                                                         ['Amodiaquine'],
                                                         ['Amodiaquine']],
                                                  dose_interval=1,
                                                  non_adherence_options=['Stop'],
                                                  non_adherence_distribution=[1],
                                                  adherence_config={
                                                        "class": "WaningEffectMapCount",
                                                        "Initial_Effect": 1,
                                                        "Durability_Map": {
                                                            "Times": [
                                                                1.0,
                                                                2.0,
                                                                3.0
                                                            ],
                                                            "Values": [
                                                                sp_resist_day1_multiply,  # for day 1
                                                                adherence,  # day 2
                                                                adherence  # day 3
                                                            ]
                                                        }
                                                    }
                                                  )
    return smc_adherent_config


def add_hfca_smc(cb, smc_df, hfca, adherence_multiplier=1, sp_resist_day1_multiply=1, seed_index=0):
    df = smc_df[smc_df['admin_name'] == hfca]
    if df.shape[0] > 0:
        if 'seed' in df.columns.values:
            df = df[df['seed'] == seed_index]
        if 'adherence' in smc_df.columns.values:
            adherent_drug_configs = smc_adherent_configuration(cb=cb,
                                                               adherence=df['adherence'].values[0]*adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        else:
            default_adherence = 0.8
            adherent_drug_configs = smc_adherent_configuration(cb=cb,
                                                               adherence=default_adherence*adherence_multiplier,
                                                               sp_resist_day1_multiply=sp_resist_day1_multiply)
        for r, row in df.iterrows() :
            cov_high_U5 = row['coverage_high_access_U5']
            cov_low_U5 = row['coverage_low_access_U5']
            cov_high_5_10 = row['coverage_high_access_5_10']
            cov_low_5_10 = row['coverage_low_access_5_10']
            if 'max_age' in smc_df.columns.values:
                max_smc_age = row['max_age']
            else:
                max_smc_age = 5
            if 'TAT' in smc_df.columns.values:
                TAT = row['TAT']
            else:
                TAT = 0

            # set diagnostic so that no one should be excluded from SMC due to fever
            add_diagnostic_survey(cb, start_day=row['simday'],
                                  coverage=1,
                                  target={"agemin": 0.25, "agemax": max(10, max_smc_age)},
                                  diagnostic_type='FEVER',
                                  diagnostic_threshold=0.5,
                                  negative_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      "class": "BroadcastEvent" }],
                                  positive_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever", # currently set so that no one should be excluded from SMC due to fever, to include fever discrimination, change to "Has_SMC_Fever"
                                      "class": "BroadcastEvent" }]
                                  )
            # U5 - high access
            add_drug_campaign(cb, 'SMC', start_days=[row['simday']],
                              coverage=cov_high_U5, target_group={'agemin': 0.25, 'agemax': 5},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'High'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U5 - low access
            add_drug_campaign(cb, 'SMC', start_days=[row['simday']],
                              coverage=cov_low_U5, target_group={'agemin': 0.25, 'agemax': 5},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'Low'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U10 (or older) - high access  (
            add_drug_campaign(cb, 'SMC', start_days=[row['simday']],
                              coverage=cov_high_5_10, target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'High'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )
            # U10 (or older) - low access
            add_drug_campaign(cb, 'SMC', start_days=[row['simday']],
                              coverage=cov_low_5_10, target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess': 'Low'}],
                              adherent_drug_configs=[adherent_drug_configs]
                              )

            if TAT:
                add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_high_U5,
                                  target_group={'agemin': 0.25, 'agemax': 5},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess' : 'High'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_low_U5,
                                  target_group={'agemin': 0.25, 'agemax': 5},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess' : 'Low'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_high_U5,
                                  target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess' : 'High'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
                add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=cov_low_U5,
                                  target_group={'agemin': 5, 'agemax': max(10, max_smc_age)},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess' : 'Low'}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
    return len(df)



def calc_high_low_access_coverages(coverage_all, high_access_frac):
    if (high_access_frac < 1) & (coverage_all >= high_access_frac):
        coverage_high = 1
        coverage_low = (coverage_all - high_access_frac) / (1 - high_access_frac)
    else:
        coverage_high = coverage_all / high_access_frac
        coverage_low = 0
    return [coverage_high, coverage_low]


def change_rtss_ips(cb):
    change_individual_property(cb,
                               target_property_name='VaccineStatus',
                               target_property_value='GotVaccine',
                               ind_property_restrictions=[{'VaccineStatus': 'None'}],
                               trigger_condition_list=['Received_Vaccine'],
                               blackout_flag=False)
    change_individual_property(cb,
                               target_property_name='VaccineStatus',
                               target_property_value='GotBooster1',
                               ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                               trigger_condition_list=['Received_Vaccine'],
                               blackout_flag=False)
    change_individual_property(cb,
                               target_property_name='VaccineStatus',
                               target_property_value='GotBooster2',
                               ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}],
                               trigger_condition_list=['Received_Vaccine'],
                               blackout_flag=False)


def add_EPI_rtss(cb, rtss_df):
    start_days = list(rtss_df['RTSS_day'].unique())
    coverage_levels = list(rtss_df['coverage'].values)
    rtss_types = list(rtss_df['vaccine'].values)
    rtss_touchpoints = list(rtss_df['rtss_touchpoints'].values)
    rtss_event_names = [f'RTSS_{x + 1}_eligible' for x in range(len(rtss_touchpoints))]

    delay_distribution_name = list(rtss_df['distribution_name'].values)[0]
    Std_Dev_list = list(rtss_df['distribution_std'].values)

    Initial_Effect_list = list(rtss_df['initial_effect'].values)
    Decay_Time_Constant_list = list(rtss_df['decay_time_constant'].values)
    try:
        Waning_Class_list = list(rtss_df['decay_class'].values)
    except:
        Waning_Class_list = ["WaningEffectExponential"] * len(rtss_touchpoints)

    for tp_time_trigger, coverage, vtype, event_name, std, init_eff, decay_t, decay_c in \
            zip(rtss_touchpoints, coverage_levels, rtss_types, rtss_event_names, Std_Dev_list,
                Initial_Effect_list, Decay_Time_Constant_list, Waning_Class_list):

        vaccine_params = {"Waning_Config": {"Initial_Effect": init_eff,
                                            "Decay_Time_Constant": decay_t,
                                            "class": decay_c}}

        if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                  "Delay_Period_Log_Normal_Mu": (
                                          np.log(tp_time_trigger + 14) - ((1 / 2) * std ** 2)),
                                  "Delay_Period_Log_Normal_Sigma": std}
            tp_time_trigger = None
        elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
            delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                  "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                  "Delay_Period_Gaussian_Std_Dev": std}
            tp_time_trigger = None
        else:
            delay_distribution = None

        # TODO: Make EPI support booster1 and booster2
        if not vtype == 'booster':
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=start_days,
                        coverage=coverage,
                        delay_distribution=delay_distribution,
                        triggered_delay=tp_time_trigger,
                        trigger_condition_list=[event_name],
                        birthtriggered=True)
        else:
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=start_days,
                        coverage=coverage,
                        delay_distribution=delay_distribution,
                        triggered_delay=tp_time_trigger,
                        trigger_condition_list=[event_name],
                        ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}],
                        birthtriggered=True)


def add_campaign_rtss(cb, rtss_df):
    for r, row in rtss_df.iterrows():
        vtype = row['vaccine']
        try:
            Waning_Class = rtss_df['rtss_decay_class_col'].unique()[0]
        except:
            Waning_Class = "WaningEffectExponential"

        vaccine_params = {"Waning_Config": {"Initial_Effect": row['initial_effect'],
                                            "Decay_Time_Constant": row['decay_time_constant'],
                                            "class": Waning_Class}}

        if vtype == 'booster1':
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=[row['RTSS_day']],
                        coverage=row['coverage'],
                        repetitions=row['repetitions'],
                        tsteps_btwn_repetitions=row['tsteps_btwn_repetitions'],
                        target_group={'agemin': row['agemin'], 'agemax': row['agemax']},
                        ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}])
        elif vtype == 'booster2':
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=[row['RTSS_day']],
                        coverage=row['coverage'],
                        repetitions=row['repetitions'],
                        tsteps_btwn_repetitions=row['tsteps_btwn_repetitions'],
                        target_group={'agemin': row['agemin'], 'agemax': row['agemax']},
                        ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}])
        elif vtype == 'booster3':
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=[row['RTSS_day']],
                        coverage=row['coverage'],
                        repetitions=row['repetitions'],
                        tsteps_btwn_repetitions=row['tsteps_btwn_repetitions'],
                        target_group={'agemin': row['agemin'], 'agemax': row['agemax']},
                        ind_property_restrictions=[{'VaccineStatus': 'GotBooster2'}])
        else:
            add_vaccine(cb,
                        vaccine_type='RTSS',
                        vaccine_params=vaccine_params,
                        start_days=[row['RTSS_day']],
                        coverage=row['coverage'],
                        repetitions=row['repetitions'],
                        tsteps_btwn_repetitions=row['tsteps_btwn_repetitions'],
                        target_group={'agemin': row['agemin'], 'agemax': row['agemax']})


def add_ds_rtss(cb, rtss_df, hfca):
    change_ips = False
    rtss_df = rtss_df[rtss_df['admin_name'].str.upper() == hfca.upper()]
    # First, process EPI style distribution
    rtss_df1 = rtss_df[rtss_df['deploy_type'] == 'EPI']
    if len(rtss_df1) > 0:
        add_EPI_rtss(cb, rtss_df1)
        change_ips = True

    # Second, process campaign style distribution
    rtss_df2 = rtss_df[rtss_df['deploy_type'] == 'campaign']
    if len(rtss_df2) > 0:
        add_campaign_rtss(cb, rtss_df2)
        change_ips = True

    if change_ips:
        change_rtss_ips(cb)

    return len(rtss_df)

def add_ds_vaccpmc(cb, pmc_df, hfca):
    df = pmc_df[pmc_df['admin_name'].str.upper() == hfca.upper()]
    if len(df) == 0:
        return 0
    try:
        num_IIV_groups = pmc_df['num_IIV_groups'].unique()[0]
    except:
        num_IIV_groups = 1

    tp_list = []
    for i, tp in enumerate(df['pmc_touchpoints']):
        tp_list.append([f'{i}', tp])
    pmc_touchpoints_dic = dict(tp_list)

    add_vaccdrug_campaign(cb, campaign_type='PMC', start_days=list(df['simday']),
                          coverages=df['coverage'],
                          target_group=pmc_touchpoints_dic,
                          delay_distribution_dic={'delay_distribution_name': df['distribution_name'],
                                                  'delay_distribution_mean': df['distribution_mean'],
                                                  'delay_distribution_std': df['distribution_std']},
                          num_IIV_groups=num_IIV_groups,
                          receiving_drugs_event=False)  ## use vaccine effects only with default offset of -10 days

    return len(pmc_df)

def add_all_interventions(cb, hfca, seed_index=1, hs_df=pd.DataFrame(), nmf_df=pd.DataFrame(), itn_df=pd.DataFrame(),
                          itn_anc_df=pd.DataFrame(), itn_use_seasonality=pd.DataFrame(), itn_decay_params=pd.DataFrame(),
                          itn_anc_adult_birthday_years=[], itn_epi_df=pd.DataFrame(),
                          itn_chw_df=pd.DataFrame(), itn_chw_annual_df=pd.DataFrame(),
                          irs_df=pd.DataFrame(), smc_df=pd.DataFrame(), pmc_df=pd.DataFrame(), vacc_df=pd.DataFrame(),
                          sp_resist_day1_multiply=1, adherence_multiplier=1, use_smc_vaccine_proxy=False):
    event_list = []

    if not irs_df.empty:
        has_irs = add_hfca_irs(cb, irs_df, hfca, seed_index=seed_index)
        if has_irs > 0:
            event_list.append('Received_IRS')
    if not smc_df.empty:
        has_smc = update_smc_access_ips(cb, smc_df=smc_df, hfca=hfca)
        if use_smc_vaccine_proxy:
            has_smc = add_hfca_vaccsmc(cb, smc_df, hfca, effective_coverage_resistance_multiplier=sp_resist_day1_multiply, seed_index=seed_index)
        else:
            has_smc = add_hfca_smc(cb, smc_df=smc_df, hfca=hfca, adherence_multiplier=adherence_multiplier, sp_resist_day1_multiply=sp_resist_day1_multiply, seed_index=seed_index)
        if has_smc > 0:
            event_list.append('Received_Campaign_Drugs')
    if not pmc_df.empty:
        has_pmc = add_ds_vaccpmc(cb, pmc_df=pmc_df, hfca=hfca)  # per default use vaccpmc
        if has_pmc > 0:
            event_list = event_list + ['Received_PMC_VaccDrug']  # 'Received_Vehicle_1','Received_Vehicle_2','Received_Vehicle_3'
    if not vacc_df.empty:
        has_vacc = add_ds_rtss(cb, rtss_df=vacc_df, hfca=hfca)
        if has_vacc > 0:
            event_list = event_list + ['Received_Vaccine']
    if not (itn_df.empty and itn_anc_df.empty and itn_epi_df.empty and itn_chw_df.empty and itn_chw_annual_df.empty):
        has_itn = add_hfca_itns(cb=cb, itn_df=itn_df, itn_anc_df=itn_anc_df,
                                itn_anc_adult_birthday_years=itn_anc_adult_birthday_years, itn_epi_df=itn_epi_df,
                                itn_chw_df=itn_chw_df, itn_chw_annual_df=itn_chw_annual_df, hfca=hfca,
                                itn_use_seasonality=itn_use_seasonality, itn_decay_params=itn_decay_params,
                                seed_index=seed_index)
        if has_itn > 0:
            event_list.append('Bednet_Got_New_One')
            event_list.append('Bednet_Using')
            event_list.append('Bednet_Discarded')
    if not hs_df.empty:
        # case management for malaria
        has_cm = add_hfca_hs(cb, hs_df, hfca, seed_index=seed_index)
        if has_cm :
            event_list.append('Received_Treatment')
            event_list.append('Received_Severe_Treatment')
        # case management for NMFs
        add_nmf_hs(cb, hs_df, nmf_df, hfca, seed_index=seed_index)
        event_list.append('Received_NMF_Treatment')
    add_event_counter_report(cb, event_trigger_list=event_list, duration=50*365)

    return {}



