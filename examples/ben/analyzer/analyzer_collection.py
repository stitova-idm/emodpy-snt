import os
import numpy as np
import pandas as pd
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
import datetime
import sys

class MonthlyPfPRAnalyzerU5(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./', start_year=2020, end_year=2023,
                 burnin=None, filter_exists=False):

        super(MonthlyPfPRAnalyzerU5, self).__init__(working_dir=working_dir,
                                                    filenames=[
                                                        f"output/MalariaSummaryReport_Monthly_{x}.json"
                                                        for x in range(start_year, end_year)]
                                                    )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.burnin = burnin
        self.filter_exists = filter_exists

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True

    def select_simulation_data(self, data, simulation):

        adf = pd.DataFrame()
        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
            pfpr = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:12]
            clinical_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:12]
            severe_cases = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:12]
            infect = [x[1] for x in d]
            d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]
            pop = [x[1] for x in d]
            d = data[fname]['DataByTime']['PfPR_2to10'][:12]
            PfPR_2to10 = d
            d = data[fname]['DataByTime']['Annual EIR'][:12]
            annualeir = d
            simdata = pd.DataFrame({'month': range(1, 13),
                                    'PfPR U5': pfpr,
                                    'Cases U5': clinical_cases,
                                    'New infections': infect,
                                    'Severe cases U5': severe_cases,
                                    'Pop U5': pop,
                                    'PfPR_2to10': PfPR_2to10,
                                    'annualeir': annualeir})
            simdata['year'] = year
            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])
            elif sweep_var == 'Run_Number' :
                adf[sweep_var] = 0

        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        print('\nSaving outputs at ' + os.path.join(self.working_dir, "U5_PfPR_ClinicalIncidence.csv"))

        adf = pd.concat(selected).reset_index(drop=True)
        if self.burnin is not None:
            adf = adf[adf['year'] > self.start_year + self.burnin]
        adf.to_csv((os.path.join(self.working_dir, 'U5_PfPR_ClinicalIncidence.csv')), index=False)

class MonthlyPfPRITNAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir="."):
        super(MonthlyPfPRITNAnalyzer, self).__init__(working_dir=working_dir,
                                                     filenames=["output/MalariaSummaryReport_Monthly_2010.json",
                                                                'output/ReportMalariaFiltered.json']
                                                     )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.mult_param = ['Habitat_Multiplier']#, 'x_Temporary_Larval_Habitat']

    # def filter(self, simulation):
    #     return simulation.tags["DS_Name_for_ITN"] == 'Dubreka'

    def select_simulation_data(self, data, simulation):

        d = data[self.filenames[0]]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
        pfpr = [x[1] for x in d]
        d = data[self.filenames[0]]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]
        age_pops = [x[1] for x in d]
        simdata = pd.DataFrame({'month': range(1, 13),
                                'PfPR U5': pfpr,
                                'Trials': age_pops})

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        all_df = pd.concat(selected).reset_index(drop=True)
        print(all_df.columns)
        grpby_var = self.mult_param + ['month', 'archetype', 'DS_Name_for_ITN']
        all_df = all_df.groupby(grpby_var)[['PfPR U5', 'Trials']].agg(
            np.mean).reset_index()
        all_df = all_df.sort_values(by=grpby_var)
        all_df.to_csv(os.path.join(self.working_dir, '%s.csv' % self.expt_name), index=False)

class MonthlyTreatedCasesAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2010, end_year=2020):
        super(MonthlyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                          filenames=["output/ReportEventCounter.json",
                                                                     "output/ReportMalariaFiltered.json"]
                                                          )
        self.sweep_variables = sweep_variables or ["LGA", "Run_Number"]
        self.channels = channels or ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

    def select_simulation_data(self, data, simulation):
        simdata = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
        simdata['Time'] = simdata.index
        if self.channels:
            d = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels})
            d['Time'] = d.index
            simdata = pd.merge(left=simdata, right=d, on='Time')
        simdata['Day'] = simdata['Time'] % 365
        simdata['Month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['Year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)
        simdata['date'] = simdata.apply(lambda x: datetime.date(int(x['Year']), int(x['Month']), 1), axis=1)

        sum_channels = ['Received_Treatment', 'Received_Severe_Treatment', 'New Clinical Cases', 'New Severe Cases',
                        'Received_NMF_Treatment']
        for x in [y for y in sum_channels if y not in simdata.columns.values]:
            simdata[x] = 0
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = simdata.groupby(['date'])[sum_channels].agg(np.sum).reset_index()
        pdf = simdata.groupby(['date'])[mean_channels].agg(np.mean).reset_index()

        simdata = pd.merge(left=pdf, right=df, on=['date'])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, 'All_Age_Monthly_Cases.csv'), index=False)
        print('\nSaving output at ' +  os.path.join(self.working_dir, 'All_Age_Monthly_Cases.csv'))


class AnnualAgebinPfPRAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./', start_year=2005,
                 end_year=2019, burnin=None):

        super(AnnualAgebinPfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                       filenames=[
                                                           f"output/MalariaSummaryReport_Annual_{start_year}to{end_year}.json"])
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.burnin = burnin

    def select_simulation_data(self, data, simulation):

        adf = pd.DataFrame()

        nyears = (self.end_year - self.start_year) + 1
        age_bins = data[self.filenames[0]]['Metadata']['Age Bins']
        pfpr2to10 = data[self.filenames[0]]['DataByTime']['PfPR_2to10'][:nyears]

        for age in range(len(age_bins)):
            d = data[self.filenames[0]]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:nyears]
            pfpr = [x[age] for x in d]
            d = data[self.filenames[0]]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:nyears]
            clinical_cases = [x[age] for x in d]
            d = data[self.filenames[0]]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:nyears]
            severe_cases = [x[age] for x in d]
            d = data[self.filenames[0]]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:nyears]
            infect = [x[age] for x in d]
            d = data[self.filenames[0]]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:nyears]
            pop = [x[age] for x in d]

            simdata = pd.DataFrame({'year': range(self.start_year, self.end_year+1),
                                    'PfPR': pfpr,
                                    'Cases': clinical_cases,
                                    'Severe cases': severe_cases,
                                    'New infections': infect,
                                    'Pop': pop})
            simdata['agebin'] = age_bins[age]
            simdata['pfpr2to10'] = pfpr2to10
            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])
            elif sweep_var == 'Run_Number' :
                adf[sweep_var] = 0

        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return
        adf = pd.concat(selected).reset_index(drop=True)

        if not os.path.exists(os.path.join(self.working_dir)):
            os.mkdir(os.path.join(self.working_dir))

        print(f'\nSaving outputs to: {self.working_dir}')

        # Discard early years used as burnin
        if self.burnin is not None:
            adf = adf[adf['year'] >= self.start_year + self.burnin]
        #adf = adf.loc[adf['agebin'] <= 100]
        adf.to_csv(os.path.join(self.working_dir, 'Agebin_PfPR_ClinicalIncidence_annual.csv'), index=False)

class MonthlyAgebinPfPRAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir='./', start_year=2020, end_year=2023,
                 burnin=None, filter_exists=False):

        super(MonthlyAgebinPfPRAnalyzer, self).__init__(working_dir=working_dir,
                                                        filenames=[
                                                            f"output/MalariaSummaryReport_Monthly_{x}.json"
                                                            for x in range(start_year, end_year)]
                                                        )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.burnin = burnin
        self.filter_exists = filter_exists

    def filter(self, simulation):
        if self.filter_exists:
            file = os.path.join(simulation.get_path(), self.filenames[0])
            return os.path.exists(file)
        else:
            return True
    def select_simulation_data(self, data, simulation):

        adf = pd.DataFrame()

        for year, fname in zip(range(self.start_year, self.end_year), self.filenames):
            age_bins = data[self.filenames[0]]['Metadata']['Age Bins']

            for age in range(len(age_bins)):
                d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]
                pfpr = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin'][:12]
                clinical_cases = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Annual Severe Incidence by Age Bin'][:12]
                severe_cases = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['New Infections by Age Bin'][:12]
                infect = [x[age] for x in d]
                d = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]
                pop = [x[age] for x in d]

                simdata = pd.DataFrame({'month': range(1, 13),
                                        'year': year,
                                        'PfPR': pfpr,
                                        'Cases': clinical_cases,
                                        'Severe cases': severe_cases,
                                        'New infections': infect,
                                        'Pop': pop})
                simdata['agebin'] = age_bins[age]
                adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                try:
                    adf[sweep_var] = simulation.tags[sweep_var]
                except:
                    adf[sweep_var] = '-'.join([str(x) for x in simulation.tags[sweep_var]])
            elif sweep_var == 'Run_Number' :
                adf[sweep_var] = 0

        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("\nWarning: No data have been returned... Exiting...")
            return
        adf = pd.concat(selected).reset_index(drop=True)

        if not os.path.exists(os.path.join(self.working_dir)):
            os.mkdir(os.path.join(self.working_dir))

        print(f'\nSaving outputs to: {self.working_dir}')

        # Discard early years used as burnin
        if self.burnin is not None:
            adf = adf[adf['year'] >= self.start_year + self.burnin]
        #adf = adf.loc[adf['agebin'] <= 100]
        adf.to_csv(os.path.join(self.working_dir, 'Agebin_PfPR_ClinicalIncidence_monthly.csv'), index=False)

class MonthlyInsetAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, rdf=pd.DataFrame(), sweep_variables=None, working_dir="."):
        super(MonthlyInsetAnalyzer, self).__init__(working_dir=working_dir,
                                                   filenames=["output/ReportMalariaFiltered.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["DS", "Run_Number"]
        self.channels = channels or ['Rainfall', 'Daily Bites per Human', 'Daily EIR', 'New Clinical Cases']
        self.rdf = rdf
        self.expt_name = expt_name

    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame( { x : data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels })
        # simdata = simdata[-365:]
        simdata['Time'] = simdata.index
        simdata['Day'] = simdata['Time'] % 365
        simdata = simdata.groupby('Day').agg(np.mean).reset_index()
        simdata['Month'] = simdata['Day'].apply(lambda x: self.monthparser((x+1) % 365))

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                if sweep_var == 'Vector_Species_Names' :
                    simdata[sweep_var] = 'funestus' if 'funestus' in simulation.tags[sweep_var] else 'no funestus'
                else :
                    simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        adf = pd.concat(selected, sort=False).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, 'monthly_inset.csv'), index=False)


class DailyBitesAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, channels=None, rdf=pd.DataFrame(), sweep_variables=None, working_dir=".",
                 firstlast = 'first'):
        super(DailyBitesAnalyzer, self).__init__(working_dir=working_dir,
                                                 filenames=["output/ReportMalariaFiltered.json"]
                                                   )
        self.sweep_variables = sweep_variables or ["DS_Name", "Run_Number"]
        self.channels = channels or ['Daily Bites per Human']
        self.rdf = rdf
        self.expt_name = expt_name
        self.firstlast = firstlast

    def filter(self, simulation):
        file = os.path.join(simulation.get_path(), self.filenames[0])
        return os.path.exists(file) 

    def select_simulation_data(self, data, simulation):
        simdata = pd.DataFrame( { x : data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels })
        if self.firstlast == 'first':
            simdata = simdata[0:365]
        else:
            simdata = simdata[-365:]
        simdata['Time'] = simdata.index

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        adf = pd.concat(selected).reset_index(drop=True)
        adf.to_csv(os.path.join(self.working_dir, 'daily_bites.csv'), index=False)


class annualSevereTreatedByAgeAnalyzer(BaseAnalyzer):

    def __init__(self, expt_name, event_name='Received_Severe_Treatment', agebins=None,
                 sweep_variables=None, working_dir=".", start_year=2010, end_year=2020,
                 ds_col='LGA'):
        super(annualSevereTreatedByAgeAnalyzer, self).__init__(working_dir=working_dir,
                                                                filenames=["output/ReportEventRecorder.csv"]
                                                                )
        self.sweep_variables = sweep_variables or ["LGA", "Run_Number"]
        self.event_name = event_name
        self.agebins = agebins or [1, 5, 125]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.ds_col = ds_col

    # def filter(self, simulation):
    #     return simulation.tags["LGA"] == 'Ushongo'

    def select_simulation_data(self, data, simulation):

        output_data = data[self.filenames[0]].copy()
        output_data = output_data[output_data['Event_Name'] == self.event_name]

        simdata = pd.DataFrame()
        if len(output_data) > 0:  # there are events of this type
            output_data['Day'] = output_data.loc[:, 'Time'] % 365
            output_data['year'] = output_data.loc[:, 'Time'].apply(lambda x: int(x / 365) + self.start_year)
            output_data['age in years'] = output_data.loc[:, 'Age'] / 365

            list_of_g = []
            for agemax in self.agebins:
                agemin = 0
                d = output_data[(output_data['age in years'] < agemax) & (output_data['age in years'] > agemin)]
                g = d.groupby(['year'])['Event_Name'].agg(len).reset_index()
                g = g.rename(columns={'Event_Name': self.event_name})
                g['agebin'] = agemax

                list_of_g = list_of_g + [g]

            simdata = pd.concat(list_of_g).reset_index(drop=True)

            for sweep_var in self.sweep_variables:
                if sweep_var in simulation.tags.keys():
                    simdata[sweep_var] = simulation.tags[sweep_var]
        else:
            simdata = pd.DataFrame(columns=['year', self.event_name] + self.sweep_variables)
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        adf = pd.concat(selected, sort=False).reset_index(drop=True)
        adf = adf.fillna(0)
        print(f"\nSaving output to {os.path.join(self.working_dir, 'Treated_Severe_Yearly_Cases_By_Age.csv')}")
        adf.to_csv(os.path.join(self.working_dir, 'Treated_Severe_Yearly_Cases_By_Age.csv'),
                   index=False)



class monthlySevereTreatedByAgeAnalyzer(BaseAnalyzer):
    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, event_name='Received_Severe_Treatment', agebins=None,
                 sweep_variables=None, working_dir=".", start_year=2010, end_year=2020,
                 ds_col='DS_Name'):
        super(monthlySevereTreatedByAgeAnalyzer, self).__init__(working_dir=working_dir,
                                                                filenames=["output/ReportEventRecorder.csv"]
                                                                )
        self.sweep_variables = sweep_variables or ["DS_Name", "Run_Number"]
        self.event_name = event_name
        self.agebins = agebins or [1, 5, 125]
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year
        self.ds_col = ds_col

    def select_simulation_data(self, data, simulation):

        output_data = data[self.filenames[0]].copy()
        output_data = output_data[output_data['Event_Name'] == self.event_name]

        simdata = pd.DataFrame()
        if len(output_data) > 0:  # there are events of this type
            output_data['Day'] = output_data.loc[:, 'Time'] % 365
            output_data['month'] = output_data.loc[:, 'Day'].apply(lambda x: self.monthparser((x + 1) % 365))
            output_data['year'] = output_data.loc[:, 'Time'].apply(lambda x: int(x / 365) + self.start_year)
            output_data['age in years'] = output_data.loc[:, 'Age'] / 365

            list_of_g = []
            for agemax in self.agebins:
                agemin = 0
                d = output_data[(output_data['age in years'] < agemax) & (output_data['age in years'] > agemin)]
                g = d.groupby(['year', 'month'])['Event_Name'].agg(len).reset_index()
                g = g.rename(columns={'Event_Name': self.event_name})
                g['agebin'] = agemax

                list_of_g = list_of_g + [g]

            simdata = pd.concat(list_of_g).reset_index(drop=True)

            for sweep_var in self.sweep_variables:
                if sweep_var in simulation.tags.keys():
                    simdata[sweep_var] = simulation.tags[sweep_var]
        else:
            simdata = pd.DataFrame(columns=['year', 'month', self.event_name] + self.sweep_variables)
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        adf = pd.concat(selected, sort=False).reset_index(drop=True)
        adf = adf.fillna(0)
        print(f"\nSaving output to {os.path.join(self.working_dir, 'Treated_Severe_Monthly_Cases_By_Age.csv')}")
        adf.to_csv(os.path.join(self.working_dir, 'Treated_Severe_Monthly_Cases_By_Age.csv'),
                   index=False)

