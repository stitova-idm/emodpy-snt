import pandas as pd
import numpy as np
from dtk.interventions.irs import add_IRS
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.property_change import change_individual_property
from malaria.interventions.adherent_drug import configure_adherent_drug
from malaria.interventions.health_seeking import add_health_seeking
from malaria.interventions.malaria_diagnostic import add_diagnostic_survey
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.interventions.malaria_vaccdrug_campaigns import add_vaccdrug_campaign
from malaria.interventions.malaria_vaccine import add_vaccine
from malaria.reports.MalariaReport import add_event_counter_report
from dtk.interventions.triggered_campaign_delay_event import triggered_campaign_delay_event

class InterventionSuite:
    # hs
    hs_ds_col = 'repDS'
    hs_duration = None 

    hs_start_col = 'simday'
    hs_coverage_age = { # column: [agemin, agemax]
        'U5_coverage': [0, 5],
        'adult_coverage': [5, 100]
    }
    hs_severe_coverage_age = {
        'severe_cases': [0, 100]
    }
    hs_rates = 0.3
    hs_severe_rates = 0.5

    # itn
    itn_ds_col = 'repDS'
    itn_leak_factor = 0.9
    itn_cov_cols = ['U5_ITN_use', 'six_nine_ITN_use', 'ten_eighteen_ITN_use', 'over_eighteen_ITN_use']
    itn_cov_age_bin = [0, 5, 10, 18]
    itn_seasonal_values = [0.032, 0.032, 0.0378, 0.154, 0.177, 0.105, 0.25, 0.32, 0.23, 0.18, 0.032]
    itn_seasonal_months = [0, 32, 60, 91, 121, 152, 182, 213, 244, 274, 364]
    itn_retention_in_yr = 1.51
    itn_preg_max_months = 48
    itn_preg_monthly_births_per = 0.00001
    itn_discard_distribution = None # 'weibull'

    # smc
    smc_ds_col = 'DS_Name'
    smc_adherence = False
    smc_default_adherence = 0.8
    smc_adherence_multiplier = 1
    smc_sp_resist_day1_multiply = 1
    smc_drug_code = 'SPA'
    smc_coverage_col = ['coverage_high_access', 'coverage_low_access'] # TODO: Consider dictionary
    smc_access = ['High', 'Low']
    smc_agemins = [0.25, 0.25]
    smc_agemax_type = ['df', 'df']
    smc_agemaxs = [5, 5]
    smc_max_age_col = 'max_age'
    smc_TAT_col = 'TAT'
    smc_leakage = True
    smc_leak_agemax = 10
    smc_leak_coverage = 0.081

    # irs
    irs_ds_col = 'DS_Name'
    irs_start_col = 'IRS_day'
    irs_coverage_col = 'coverage'
    irs_box_dur_col = 'box_duration'
    irs_decay_t_col = 'decay_time_constant'
    irs_init_eff_col = 'initial_killing'

    # rtss
    rtss_ds_col = 'DS_Name'
    rtss_type_col = 'vaccine'
    rtss_start_col = 'RTSS_day'
    rtss_coverage_col = 'coverage_levels'
    rtss_touchpoint_col = 'rtss_touchpoints'  # days since births!
    rtss_distribution_col = 'distribution_name'
    rtss_std_col = 'distribution_std'
    rtss_min_age_col = 'agemin'
    rtss_max_age_col = 'agemax'
    rtss_repetitions = 'repetitions'
    rtss_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'
    rtss_init_eff_col = 'initial_killing'
    rtss_decay_t_col = 'decay_time_constant'
    rtss_decay_class_col = 'decay_class'
    rtss_auto_changeips = True

    # PMC
    pmc_ds_col = 'DS_Name'
    pmc_start_col = 'PMC_day'
    pmc_correlated_cov_col = 'pmc_correlated_cov'  # boolean, requires PMC IP in demographics file
    pmc_coverage_col = 'coverage_levels'
    pmc_touchpoint_col = 'pmc_touchpoints'  # days since births!
    pmc_agemin = 'agemin'
    pmc_agemax = 'agemax'
    pmc_repetitions = 'repetitions'
    pmc_tsteps_btwn_repetitions = 'tsteps_btwn_repetitions'

    def add_ds_hs(self, cb, hs_df, my_ds) :
        ds_col = self.hs_ds_col
        duration = self.hs_duration 
        df = hs_df[hs_df[ds_col] == my_ds]
        for r, row in df.iterrows() :
            self.add_hs_from_file(cb, row, duration=duration) 

        return len(df)


    def add_hs_from_file(self, cb, row, duration): 
        rates = self.hs_rates
        severe_rates = self.hs_severe_rates

        start_day = row[self.hs_start_col] #if start_day_override < 0 else start_day_override
        if duration is None:
            duration = row['duration']
        
        # Age-specific case management coverage, per default if custom_age_coverage in hs_df
        if 'custom_age_coverage' in row.index:
            self.hs_coverage_age= {'custom_age_coverage': [int(row['agemin']), int(row['agemax'])]}
            self.hs_severe_coverage_age = {'severe_cases': [int(row['agemin']), int(row['agemax'])]}

        # Uncomplicated
        targets = []
        for key, value in self.hs_coverage_age.items():
            targets.append({
                               'trigger': 'NewClinicalCase',
                               'coverage': row[key],
                               'agemin': value[0],
                               'agemax': value[1],
                               'seek': 1,
                               'rate': rates
                           }) 
        add_health_seeking(cb, start_day=start_day,
                           targets=targets,
                           drug=['Artemether', 'Lumefantrine'], duration=duration)
        
        # Severe
        targets = []
        for key, value in self.hs_severe_coverage_age.items():
            targets.append({
                               'trigger': 'NewSevereCase',
                               'coverage': row[key],
                               'agemin': value[0],
                               'agemax': value[1],
                               'seek': 1,
                               'rate': severe_rates
                           }) 
        add_health_seeking(cb, start_day=start_day,
                           targets=targets,
                           drug=['Artemether', 'Lumefantrine'], duration=duration,
                           broadcast_event_name='Received_Severe_Treatment')
    
    def add_ds_itns(self, cb, itn_df, my_ds) :
        itn_ds_col = self.itn_ds_col
        df = itn_df[itn_df[itn_ds_col].str.upper() == my_ds.upper()]
        df = df.drop_duplicates()
        nets = len(df)
        for r, row in df.iterrows() :
            if 'blocking_rate' not in row.index:
                row['blocking_rate'] = row['block_initial'] #FIXME in input files for nigeria
            self.add_itn_by_row(cb, row)

        return nets

    def add_ds_itns_addtnl(self, cb, itn_addtnl_df, my_ds):
        itn_ds_col = self.itn_ds_col

        df = itn_addtnl_df[itn_addtnl_df[itn_ds_col].str.upper() == my_ds.upper()]
        df = df.drop_duplicates()
        nets = len(df)
        for r, row in df.iterrows() :
            itn_type = row['type']
            if 'blocking_rate' not in row.index:
                row['blocking_rate'] = row['block_initial'] #FIXME in input files for nigeria
            func = getattr(self, 'add_itn_' + itn_type + '_by_row')
            func(cb, row)
        
        return nets

    def add_itn_by_row(self, cb, row) :
        usages = [row[x] for x in self.itn_cov_cols]
        coverage_all = np.max(usages)
        if coverage_all == 0:
            coverage_all = 1
        usages = [x/coverage_all for x in usages]
        
        self.ITN_age_season_template(cb, start=row['simday'],
                                     demographic_coverage=coverage_all,
                                     killing_rate=row['kill_rate'],
                                     blocking_rate=row['blocking_rate'],
                                     age_bin=self.itn_cov_age_bin,
                                     usages=[x * self.itn_leak_factor for x in usages])

    def ITN_age_season_template(self, cb, start, demographic_coverage,
                                killing_rate, blocking_rate,
                                age_bin, usages, duration=-1, 
                                birth_triggered=False,
                                trigger_condition_list=None,
                                ind_property_restrictions=None):
        seasonal_times, seasonal_scales = self.adjust_itn_seasonals(start)

        if self.itn_discard_distribution == 'weibull':
            discard_config = {"Expiration_Period_Distribution": "WEIBULL_DISTRIBUTION",
                              "Expiration_Period_Kappa": 2.3,
                              "Expiration_Period_Lambda": 2.5 * 365}
        else:
            discard_config = {"Expiration_Period_Distribution": "DUAL_EXPONENTIAL_DISTRIBUTION",
                                           "Expiration_Period_Proportion_1": 0.9,
                                           "Expiration_Period_Mean_1": 365*self.itn_retention_in_yr, # Burkina 1.7
                                           "Expiration_Period_Mean_2": 3650}

        add_ITN_age_season(cb, start=start,
                           demographic_coverage=demographic_coverage,
                           killing_config={
                               "Initial_Effect": killing_rate,
                               "Decay_Time_Constant": 1460,
                               "class": "WaningEffectExponential"},
                           blocking_config={
                               "Initial_Effect": blocking_rate,
                               "Decay_Time_Constant": 730,
                               "class": "WaningEffectExponential"},
                           discard_times=discard_config,
                           age_dependence={'Times': age_bin,
                                           'Values': usages},
                           seasonal_dependence={"Times": seasonal_times, "Values": seasonal_scales},
                           duration=duration, birth_triggered=birth_triggered,
                           trigger_condition_list=trigger_condition_list,
                           ind_property_restrictions=ind_property_restrictions
                           )

    def adjust_itn_seasonals(self, simday):
        seasonal_scales = [x/max(self.itn_seasonal_values) for x in self.itn_seasonal_values]

        seasonal_offset = simday % 365
        seasonal_times = [(x + (365 - seasonal_offset)) % 365 for x in self.itn_seasonal_months]

        zipped_lists = zip(seasonal_times, seasonal_scales)
        sorted_pairs = sorted(zipped_lists)
        tuples = zip(*sorted_pairs)
        seasonal_times, seasonal_scales = [list(tuple) for tuple in tuples]
        if seasonal_times[0] > 0 :
            seasonal_times.insert(0, 0)
            seasonal_scales.insert(0, seasonal_scales[-1])

        return(seasonal_times, seasonal_scales)
    
    def add_itn_antenatal_by_row(self, cb, row) :
        # for r, row in itn_anc_df.iterrows() :
        self.ITN_age_season_template(cb, start=row['simday'],
                                     demographic_coverage=row['coverage'],
                                     killing_rate=row['kill_rate'],
                                     blocking_rate=row['blocking_rate'],
                                     age_bin=[0, 5],
                                     usages=[self.itn_leak_factor for i in range(2)],
                                     birth_triggered=True)

    def add_itn_pregnant_by_row(self, cb, row) :
        for mm in range(self.itn_preg_max_months):
            self.ITN_age_season_template(cb, start=30*mm,
                                         demographic_coverage=row['coverage'],
                                         killing_rate=row['kill_rate'],
                                         blocking_rate=row['blocking_rate'],
                                         age_bin=[0, 5],
                                         usages=[self.itn_leak_factor for i in range(2)],
                                         ind_property_restrictions=[{'AgeGroup': '15to30'}])

    def add_itn_severe_by_row(self, cb, row) :
        usages = [row[x] for x in self.itn_cov_cols]
        coverage_all = np.max(usages)
        if coverage_all == 0:
            coverage_all = 1
        usages = [x/coverage_all for x in usages]
        
        self.ITN_age_season_template(cb, start=0,
                                     demographic_coverage=row['coverage'],
                                     killing_rate=row['kill_rate'],
                                     blocking_rate=row['blocking_rate'],
                                     age_bin=self.itn_cov_age_bin,
                                     usages=[x * self.itn_leak_factor for x in usages],
                                     trigger_condition_list=['Received_Severe_Treatment'])
   
    def add_itn_uncomplicated_by_row(self, cb, row) :
        usages = [row[x] for x in self.itn_cov_cols]
        coverage_all = np.max(usages)
        if coverage_all == 0:
            coverage_all = 1
        usages = [x/coverage_all for x in usages]
        
        self.ITN_age_season_template(cb, start=0,
                                     demographic_coverage=row['coverage'],
                                     killing_rate=row['kill_rate'],
                                     blocking_rate=row['blocking_rate'],
                                     age_bin=self.itn_cov_age_bin,
                                     usages=[x * self.itn_leak_factor for x in usages],
                                     trigger_condition_list=['Received_Treatment'])

    def add_ds_vaccsmc(self, cb, smc_df, my_ds):
        ds_col = self.smc_ds_col
        # agemax_type = 'fixed' determined by agemaxs, 'df' determined by row['max_age']
        df = smc_df[smc_df[ds_col] == my_ds].reset_index()
        if len(df) == 0: 
            return len(df)

        def det_agemax(type, dfmax, argmax):
            if type == 'fixed':
                return (argmax)
            elif type == 'df':
                return (dfmax)

        if self.smc_max_age_col in df.columns.values:
            max_smc_age = list(df[self.smc_max_age_col])[0]  # assume same max age for all rounds
        else:
            max_smc_age = 5
        if self.smc_TAT_col in df.columns.values:
            TAT = list(df[self.smc_TAT_col])[0]  # assume same for all rounds
        else:
            TAT = 0  # assume kids with fever don't get SMC at all and are also not referred

        coverage_col = self.smc_coverage_col
        agemins = self.smc_agemins
        agemaxs = self.smc_agemaxs
        agemax_type = self.smc_agemax_type
        access = self.smc_access
        agemax_actual = []
        for i in range(len(coverage_col)):
            agemax_actual.append(det_agemax(agemax_type[i], max_smc_age, agemaxs[i]))
            add_vaccdrug_campaign(cb, campaign_type='SMC', start_days=df['simday'],
                                  coverages=df[coverage_col[i]],
                                  target_group={'agemin': agemins[i],
                                                'agemax': agemax_actual[i]},
                                  #trigger_condition_list=['No_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': access[i]}],
                                  receiving_drugs_event = False)  ## If False uses vaccSMC with automatic offset of 17 days, if True, uses vaccDrugSMC
                                  
            if TAT:  # TODO
                raise ValueError('TAT not yet implemented for vaccSMC')

        #for r, row in df.iterrows():
        #    add_diagnostic_survey(cb, start_day=row['simday'],
        #                          coverage=1,
        #                          target={"agemin": np.min(agemins),
        #                                  "agemax": np.max(agemax_actual)},
        #                          diagnostic_type='FEVER',
        #                          diagnostic_threshold=0.5,
        #                          negative_diagnosis_configs=[{
        #                              "Broadcast_Event": "No_SMC_Fever",
        #                              "class": "BroadcastEvent"}],
        #                          positive_diagnosis_configs=[{
        #                              "Broadcast_Event": "Has_SMC_Fever",
        #                              "class": "BroadcastEvent"}]
        #                          )

        if self.smc_leakage and len(df) > 0:
            # leakage to between age max of SMC group and user defined max
            leak_agemin = np.max(agemax_actual)
            smc_leak_coverages = [self.smc_leak_coverage] * len(df['simday'])
            add_vaccdrug_campaign(cb, campaign_type='SMC', start_days=df['simday'],
                                  coverages=smc_leak_coverages,
                                  target_group={'agemin': leak_agemin,
                                                'agemax': self.smc_leak_agemax},
                                  #trigger_condition_list=['No_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess': access[i]}])

        return len(df)


    def add_ds_smc(self, cb, smc_df, my_ds) :
        ds_col = self.smc_ds_col
        adherence_multiplier = self.smc_adherence_multiplier
        sp_resist_day1_multiply = self.smc_sp_resist_day1_multiply
        # agemax_type = 'fixed' determined by agemaxs, 'df' determined by row['max_age']
        df = smc_df[smc_df[ds_col] == my_ds]
        drug_code = self.smc_drug_code
        if self.smc_adherence:
            if 'adherence' in smc_df.columns.values:
                adherent_drug_configs = self.smc_adherent_configuration(cb=cb,
                                                                        adherence=df['adherence'].values[0]*adherence_multiplier,
                                                                        sp_resist_day1_multiply=sp_resist_day1_multiply)
            else:
                default_adherence = self.smc_default_adherence
                adherent_drug_configs = self.smc_adherent_configuration(cb=cb,
                                                                        adherence=default_adherence*adherence_multiplier,
                                                                        sp_resist_day1_multiply=sp_resist_day1_multiply)
            drug_code = None
            adherent_drug_configs = [adherent_drug_configs]

        else:
            adherent_drug_configs = None

        def det_agemax (type, dfmax, argmax):
            if type == 'fixed':
                return(argmax)
            elif type == 'df':
                return(dfmax)

        for r, row in df.iterrows() : 
            if self.smc_max_age_col in smc_df.columns.values:
                max_smc_age = row[self.smc_max_age_col]
            else:
                max_smc_age = 5
            if self.smc_TAT_col in smc_df.columns.values:
                TAT = row[self.smc_TAT_col]
            else:
                TAT = 0 # assume kids with fever don't get SMC at all and are also not referred

            coverage_col = self.smc_coverage_col
            agemins = self.smc_agemins
            agemaxs = self.smc_agemaxs
            agemax_type = self.smc_agemax_type
            access = self.smc_access
            agemax_actual = []
            for i in range(len(coverage_col)):
                agemax_actual.append(det_agemax(agemax_type[i], max_smc_age, agemaxs[i]))
                add_drug_campaign(cb, 'SMC', drug_code, start_days=[row['simday']],
                              coverage=row[coverage_col[i]], 
                              target_group={'agemin' : agemins[i], 
                                            'agemax' : agemax_actual[i]},
                              listening_duration=2,
                              trigger_condition_list=['No_SMC_Fever'],
                              ind_property_restrictions=[{'SMCAccess' : access[i]}],
                              adherent_drug_configs=adherent_drug_configs)
                if TAT:
                    add_drug_campaign(cb, 'MDA', drug_code='AL', start_days=[row['simday']],
                                  coverage=row[coverage_col[i]],
                                  target_group={'agemin': agemins[i], 
                                                'agemax': agemax_actual[i]},
                                  listening_duration=2,
                                  trigger_condition_list=['Has_SMC_Fever'],
                                  ind_property_restrictions=[{'SMCAccess' : access[i]}],
                                  receiving_drugs_event_name='Received_TAT_Treatment')
        
            add_diagnostic_survey(cb, start_day=row['simday'],
                                  coverage=1,
                                  target={"agemin": np.min(agemins), 
                                          "agemax": np.max(agemax_actual)},
                                  diagnostic_type='FEVER',
                                  diagnostic_threshold=0.5,
                                  negative_diagnosis_configs=[{
                                      "Broadcast_Event": "No_SMC_Fever",
                                      "class": "BroadcastEvent" }],
                                  positive_diagnosis_configs=[{
                                      "Broadcast_Event": "Has_SMC_Fever",
                                      "class": "BroadcastEvent" }]
                                  )

        if self.smc_leakage and len(df) > 0:
            # leakage to between age max of SMC group and user defined max
            leak_agemin = np.max(agemax_actual)
            add_drug_campaign(cb, 'SMC', drug_code, start_days=df['simday'].values, # Might throw error?
                              coverage=self.smc_leak_coverage, 
                              target_group={'agemin' : leak_agemin, 'agemax' : self.smc_leak_agemax})

        return len(df)

    def smc_adherent_configuration(self, cb, adherence, sp_resist_day1_multiply):
        smc_adherent_config = configure_adherent_drug(cb, 
                                                      doses = [["Sulfadoxine", "Pyrimethamine",'Amodiaquine'],
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

    def add_ds_irs(self, cb, irs_df, my_ds) :
        irs_df = irs_df[irs_df[self.irs_ds_col].str.upper() == my_ds.upper()]
        for r, row in irs_df.iterrows() :
            add_IRS(cb, start=row[self.irs_start_col], 
                    coverage_by_ages=[{"coverage": row[self.irs_coverage_col], "min": 0, "max": 100}],
                    killing_config={
                        "class": "WaningEffectBoxExponential",
                        "Box_Duration": row[self.irs_box_dur_col],  # based on PMI data from Burkina
                        "Decay_Time_Constant": row[self.irs_decay_t_col],
                        "Initial_Effect": row[self.irs_init_eff_col]},
                    )

        return len(irs_df)

    def change_rtss_ips (self, cb):
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


    def add_EPI_rtss(self, cb, rtss_df):
        start_days = list(rtss_df[self.rtss_start_col].unique())
        coverage_levels = list(rtss_df[self.rtss_coverage_col].values)
        rtss_types = list(rtss_df[self.rtss_type_col].values)
        rtss_touchpoints = list(rtss_df[self.rtss_touchpoint_col].values)
        rtss_event_names = [f'RTSS_{x + 1}_eligible' for x in range(len(rtss_touchpoints))]
        
        delay_distribution_name = list(rtss_df[self.rtss_distribution_col].values)[0]
        Std_Dev_list = list(rtss_df[self.rtss_std_col].values)
        
        Initial_Effect_list = list(rtss_df[self.rtss_init_eff_col].values)
        Decay_Time_Constant_list = list(rtss_df[self.rtss_decay_t_col].values)
        try:
            Waning_Class_list = list(rtss_df[self.rtss_decay_class_col].values)
        except:
            Waning_Class_list = ["WaningEffectExponential"] * len(rtss_touchpoints)
        
        for tp_time_trigger, coverage, vtype, event_name, std, init_eff,decay_t,decay_c in \
                zip(rtss_touchpoints, coverage_levels, rtss_types, rtss_event_names,Std_Dev_list,
                    Initial_Effect_list,Decay_Time_Constant_list,Waning_Class_list):
        
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

    def add_campaign_rtss(self, cb, rtss_df):
        for r, row in rtss_df.iterrows():
            vtype = row[self.rtss_type_col]
            try:
                Waning_Class = rtss_df['rtss_decay_class_col'].unique()[0]
            except:
                Waning_Class = "WaningEffectExponential"
     
            vaccine_params = {"Waning_Config": {"Initial_Effect": row[self.rtss_init_eff_col],
                                                "Decay_Time_Constant": row[self.rtss_decay_t_col],
                                                "class": Waning_Class}}
     
            if vtype == 'booster1':
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=[row[self.rtss_start_col]],
                            coverage=row[self.rtss_coverage_col],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'VaccineStatus': 'GotVaccine'}])
            elif vtype == 'booster2':
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=[row[self.rtss_start_col]],
                            coverage=row[self.rtss_coverage_col],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'VaccineStatus': 'GotBooster1'}])
            elif vtype == 'booster3':
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=[row[self.rtss_start_col]],
                            coverage=row[self.rtss_coverage_col],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]},
                            ind_property_restrictions=[{'VaccineStatus': 'GotBooster2'}])
            else:
                add_vaccine(cb,
                            vaccine_type='RTSS',
                            vaccine_params=vaccine_params,
                            start_days=[row[self.rtss_start_col]],
                            coverage=row[self.rtss_coverage_col],
                            repetitions=row[self.rtss_repetitions],
                            tsteps_btwn_repetitions=row[self.rtss_tsteps_btwn_repetitions],
                            target_group={'agemin': row[self.rtss_min_age_col], 'agemax': row[self.rtss_max_age_col]})


    def add_ds_rtss(self, cb, rtss_df, my_ds):
        rtss_df = rtss_df[rtss_df[self.rtss_ds_col].str.upper() == my_ds.upper()]
        # First, process EPI style distribution
        rtss_df1 = rtss_df[rtss_df['deploy_type'] == 'EPI']
        if len(rtss_df1) > 0:
            self.add_EPI_rtss(cb, rtss_df1)
        
        # Second, process campaign style distribution
        rtss_df2 = rtss_df[rtss_df['deploy_type'] == 'campaign']
        if len(rtss_df2) > 0:
            self.add_campaign_rtss(cb, rtss_df2)

        if self.rtss_auto_changeips:
            self.change_rtss_ips(cb)

        return len(rtss_df)

    def add_ds_vaccpmc(self, cb, pmc_df, my_ds):
        df = pmc_df[pmc_df[self.pmc_ds_col].str.upper() == my_ds.upper()]
        if len(df) == 0: return 0
        try:
            num_IIV_groups= pmc_df['num_IIV_groups'].unique()[0]
        except:
            num_IIV_groups=1

        tp_list = []
        for i, tp in enumerate(df[self.pmc_touchpoint_col]):
            tp_list.append([f'{i}', tp])
        pmc_touchpoints_dic = dict(tp_list)

        add_vaccdrug_campaign(cb, campaign_type='PMC', start_days=list(df[self.pmc_start_col]),
                              coverages=df[self.pmc_coverage_col],
                              target_group=pmc_touchpoints_dic,
                              delay_distribution_dic={'delay_distribution_name': df['distribution_name'],
                                                      'delay_distribution_mean': df['distribution_mean'],
                                                      'delay_distribution_std': df['distribution_std']},
                              num_IIV_groups=num_IIV_groups,
                              receiving_drugs_event = False) ## use vaccine effects only with default offset of -10 days

        return len(pmc_df)

    def add_ds_pmc(self, cb, pmc_df, my_ds):
        pmc_df = pmc_df[pmc_df[self.pmc_ds_col].str.upper() == my_ds.upper()]
        EPI = pmc_df['deploy_type'].unique()[0] == 'EPI'
        drug_code = pmc_df['drug_code'].unique()[0]

        try:
            num_IIV_groups= pmc_df['num_IIV_groups'].unique()[0]
            pmc_IIV= True
        except:
            num_IIV_groups=None
            pmc_IIV = False

        if not EPI:
            """Use campaign-style deployment"""
            for r, row in pmc_df.iterrows():
                add_drug_campaign(cb, campaign_type='SMC',
                                  drug_code=drug_code,
                                  coverage=row[self.pmc_coverage_col],
                                  start_days=[row[self.pmc_start_col]],
                                  target_group={'agemin': row[self.pmc_agemin], 'agemax': row[self.pmc_agemax]},
                                  repetitions=row[self.pmc_repetitions],
                                  tsteps_btwn_repetitions=row[self.pmc_tsteps_btwn_repetitions],
                                  receiving_drugs_event_name='Received_PMC')

        else:
            """Use birth-triggered deployment for PMC administered along EPI"""
            try:
                pmc_correlated_cov = bool(pmc_df[self.pmc_correlated_cov_col].values[0])
            except:
                pmc_correlated_cov = False

            coverage_levels = list(pmc_df[self.pmc_coverage_col].values)
            pmc_touchpoints = list(pmc_df[self.pmc_touchpoint_col].values)
            start_days = list(pmc_df[self.pmc_start_col].unique())
            delay_distribution_name = list(pmc_df['distribution_name'].values)[0]
            Std_Dev_list = list(pmc_df['distribution_std'].values)
            pmc_event_names = [f'PMC_{x + 1}' for x in range(len(pmc_touchpoints))]

            if not pmc_IIV:
                for tp_time_trigger, coverage, event_name, std in zip(pmc_touchpoints, coverage_levels,
                                                                      pmc_event_names,Std_Dev_list):

                    if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
                        delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Delay_Period_Log_Normal_Mu": (np.log(tp_time_trigger + 14) - ((1 / 2) * std ** 2)),
                                              "Delay_Period_Log_Normal_Sigma": std}
                        tp_time_trigger = None
                    elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
                        delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                              "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                              "Delay_Period_Gaussian_Std_Dev": std}
                        tp_time_trigger = None
                    else:
                        delay_distribution = None


                    if pmc_correlated_cov:
                        add_drug_campaign(cb,
                                          campaign_type='PMC',
                                          drug_code=drug_code,
                                          start_days=start_days,
                                          coverage=coverage, # align with IP PMC subgroup in demographics file
                                          delay_distribution=delay_distribution,
                                          triggered_campaign_delay=tp_time_trigger,
                                          trigger_name=event_name,
                                          ind_property_restrictions=[{'PMCEligibility':'PMC'}])
                    else:
                        add_drug_campaign(cb,
                                          campaign_type='PMC',
                                          drug_code=drug_code,
                                          start_days=start_days,
                                          coverage=coverage,
                                          delay_distribution=delay_distribution,
                                          triggered_campaign_delay=tp_time_trigger,
                                          trigger_name=event_name)
            else:
                """When running with IIV birth triggered event needs to happen separate from drug campaigns"""
                IIV_groups = ["Group%d" % x for x in range(num_IIV_groups)]
                for tp_time_trigger, coverage, event_name, std in zip(pmc_touchpoints, coverage_levels,
                                                                      pmc_event_names,Std_Dev_list):

                    if delay_distribution_name == "LOG_NORMAL_DISTRIBUTION":
                        delay_distribution = {"Delay_Period_Distribution": "LOG_NORMAL_DISTRIBUTION",
                                              "Delay_Period_Log_Normal_Mu": (np.log(tp_time_trigger + 14) - ((1 / 2) * std ** 2)),
                                              "Delay_Period_Log_Normal_Sigma": std}
                        tp_time_trigger = None
                    elif delay_distribution_name == "GAUSSIAN_DISTRIBUTION":
                        delay_distribution = {"Delay_Period_Distribution": "GAUSSIAN_DISTRIBUTION",
                                              "Delay_Period_Gaussian_Mean": tp_time_trigger,
                                              "Delay_Period_Gaussian_Std_Dev": std}
                        tp_time_trigger = None
                    else:
                        delay_distribution = None

                    triggered_campaign_delay_event(cb, start=start_days[0],
                                                   delay_distribution=delay_distribution,
                                                   triggered_campaign_delay=tp_time_trigger,
                                                   coverage=1,
                                                   trigger_condition_list=['Births'],
                                                   event_to_send_out=event_name)

                    for index, val in enumerate(IIV_groups):
                        if drug_code == 'SDX_PYR':
                            drug_doses = [['SulfadoxinePyrimethamine_%d' % index]]
                        elif drug_code == 'SP':
                            drug_doses = [['Sulfadoxine_%d' % index,
                                           'Pyrimethamine_%d' % index]]
                        else:
                            raise ValueError('drug_code for PMC needs to be "SP" or "SDX_PYR"')

                        dose_response_config = configure_adherent_drug(cb,
                                                                       doses=drug_doses,
                                                                       dose_interval=1,
                                                                       non_adherence_options=['Stop'],
                                                                       non_adherence_distribution=[1],
                                                                       adherence_config={
                                                                           "class": "WaningEffectMapCount",
                                                                           "Initial_Effect": 1,
                                                                           "Durability_Map": {
                                                                               "Times": [
                                                                                   1.0
                                                                               ],
                                                                               "Values": [
                                                                                   1
                                                                               ]
                                                                           }
                                                                       }
                                                                       )

                        add_drug_campaign(cb, campaign_type='MDA',
                                          adherent_drug_configs=[dose_response_config],
                                          coverage=coverage, start_days=start_days,
                                          ind_property_restrictions=[{'DrugResponseGroup': val}],
                                          repetitions=1,
                                          delay_distribution=delay_distribution,
                                          receiving_drugs_event_name=f'Received_{event_name}',
                                          trigger_condition_list=[event_name])

        return len(pmc_df)

def add_all_interventions(cb, int_suite, my_ds, hs_df=pd.DataFrame(), 
                          itn_df=pd.DataFrame(),  
                          itn_addtnl_df=pd.DataFrame(), 
                          smc_df=pd.DataFrame(),
                          irs_df=pd.DataFrame(), 
			              rtss_df=pd.DataFrame(),
                          pmc_df=pd.DataFrame(),
                          addtl_smc_func=None,
                          nmf_as_default=True,
                          add_event_report=True,
                          rep_start=0,
                          rep_duration=10000): # itn_anc_df, irs_df
    if nmf_as_default:
        event_list = ['Received_NMF_Treatment']
    else:
        event_list = []
    if not irs_df.empty:
       has_irs = int_suite.add_ds_irs(cb, irs_df, my_ds)
       if has_irs > 0:
            event_list.append('Received_IRS')

    if not smc_df[smc_df[int_suite.smc_ds_col] == my_ds].empty :
        if addtl_smc_func:
            addtl_smc = addtl_smc_func(cb, smc_df, my_ds)
        # has_smc = int_suite.add_ds_smc(cb, smc_df, my_ds)  #  'Received_Campaign_Drugs'
        has_smc = int_suite.add_ds_vaccsmc(cb, smc_df, my_ds)  # per default use vaccsmc
        if has_smc > 0:
            event_list = event_list + ['Received_SMC_VaccDrug'] # 'Received_Vehicle'

    if not itn_df.empty:
        has_itn = int_suite.add_ds_itns(cb, itn_df, my_ds)
        if has_itn > 0:
            event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']

    if not itn_addtnl_df.empty:
        has_itn_addtnl = int_suite.add_ds_itns_addtnl(cb, itn_addtnl_df, my_ds)
        if has_itn_addtnl > 0: 
            event_list = event_list + ['Bednet_Got_New_One', 'Bednet_Using', 'Bednet_Discarded']

    if not hs_df.empty :
        has_cm = int_suite.add_ds_hs(cb, hs_df, my_ds)
        if has_cm :
            event_list.append('Received_Treatment')
            event_list.append('Received_Severe_Treatment')
    
    if not rtss_df.empty:
        has_rtss = int_suite.add_ds_rtss(cb, rtss_df, my_ds)
        if has_rtss > 0:
            event_list.append('Received_Vaccine')											 

    if not pmc_df.empty:
        # has_pmc = int_suite.add_ds_pmc(cb, pmc_df, my_ds)  #  'Received_PMC'
        has_pmc = int_suite.add_ds_vaccpmc(cb, pmc_df, my_ds) # per default use vaccpmc
        if has_pmc > 0:
            event_list = event_list + ['Received_PMC_VaccDrug']  # 'Received_Vehicle_1','Received_Vehicle_2','Received_Vehicle_3'

    event_list = list(np.unique(event_list))
    if len(event_list) > 0 and add_event_report:
        add_event_counter_report(cb, event_trigger_list=event_list, start=rep_start,
                                 duration=rep_duration)
    return {}


def update_smc_access_ips(cb, smc_df, my_ds) :
    df = smc_df[smc_df['DS_Name'] == my_ds]

    # copied directly from burkina set_up_simulation_config
    # change SMCAccess property of newborns (important for those born after the first SMC round in a year)
    change_individual_property(cb, 'SMCAccess', 'Low', target_group={'agemin': 0, 'agemax': 5}, coverage=1, blackout_flag=False)
    change_individual_property(cb, 'SMCAccess', 'High', target_group={'agemin': 0, 'agemax': 5}, coverage=df['high_access_U5'].values[0], blackout_flag=False)

    # before the first SMC round in each year, change the SMCAccess IP for the U5 and O5 age groups
    # simdays of the first rounds (change IPs one week before)
    first_round_days = df.loc[df['round'] == 1, 'simday'].values
    change_IP_days = [first_round_days[yy] - 7 for yy in range(len(first_round_days))]

    for rr in change_IP_days:
        change_individual_property(cb, 'SMCAccess', 'Low', target_group={'agemin': 0, 'agemax': 5}, start_day=rr, coverage=1, blackout_flag=False)
        change_individual_property(cb, 'SMCAccess', 'High', target_group={'agemin': 0, 'agemax': 5}, start_day=rr, coverage=df['high_access_U5'].values[0], blackout_flag=False)
        change_individual_property(cb, 'SMCAccess', 'Low', target_group={'agemin': 5, 'agemax': 120}, start_day=rr, coverage=1, blackout_flag=False)
        change_individual_property(cb, 'SMCAccess', 'High', target_group={'agemin': 5, 'agemax': 120}, start_day=rr, coverage=df['high_access_5_10'].values[0], blackout_flag=False)

    return {'DS_Name' : my_ds}
