import pandas as pd
import numpy as np
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
import datetime
import os
import sys
sys.path.append('../')



class monthlyTreatedCasesAnalyzer(BaseAnalyzer):

    @classmethod
    def monthparser(self, x):
        if x == 0:
            return 12
        else:
            return datetime.datetime.strptime(str(x), '%j').month

    def __init__(self, expt_name, channels=None, sweep_variables=None, working_dir=".", start_year=2020, end_year=2026):
        super(monthlyTreatedCasesAnalyzer, self).__init__(working_dir=working_dir,
                                                          filenames=["output/ReportEventCounter.json",
                                                                     "output/ReportMalariaFiltered.json"]
                                                          )
        self.sweep_variables = sweep_variables or ["__sample_index__", "Run_Number"]
        if channels is None:
            self.channels = ['Received_Treatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment']
        else:
            self.channels = channels
        self.inset_channels = ['Statistical Population', 'New Clinical Cases', 'New Severe Cases', 'PfHRP2 Prevalence']
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year

        # added to bypass failed cases
        # def filter(self, simulation):
        #     return simulation.status.name == 'Succeeded'

    def filter(self, simulation):
        # return simulation.tags['__sample_index__'] == 0  #!!!
        return simulation

    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame({x: data[self.filenames[0]]['Channels'][x]['Data'] for x in self.channels})
        simdata['Time'] = simdata.index

        d = pd.DataFrame({x: data[self.filenames[1]]['Channels'][x]['Data'] for x in self.inset_channels})
        d['Time'] = d.index

        if len(self.channels) > 0:
            simdata = pd.merge(left=simdata, right=d, on='Time')
        else:
            simdata = d
        simdata['Day'] = simdata['Time'] % 365
        simdata['month'] = simdata['Day'].apply(lambda x: self.monthparser((x + 1) % 365))
        simdata['year'] = simdata['Time'].apply(lambda x: int(x / 365) + self.start_year)

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                simdata[sweep_var] = simulation.tags[sweep_var]
        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(os.path.join(self.working_dir, self.expt_name)):
            os.mkdir(os.path.join(self.working_dir, self.expt_name))

        adf = pd.concat(selected).reset_index(drop=True)
        adf['date'] = adf.apply(lambda x: datetime.date(int(x['year']), int(x['month']), 1), axis=1)

        sum_channels = self.channels + ['New Clinical Cases', 'New Severe Cases']
        mean_channels = ['Statistical Population', 'PfHRP2 Prevalence']

        df = adf.groupby(['__sample_index__', 'date', 'Run_Number'])[sum_channels].agg(np.sum).reset_index()
        pdf = adf.groupby(['__sample_index__', 'date', 'Run_Number'])[mean_channels].agg(np.mean).reset_index()

        adf = pd.merge(left=pdf, right=df, on=['__sample_index__', 'date', 'Run_Number'])
        adf.to_csv(os.path.join(self.working_dir, self.expt_name, 'All_Age_monthly_Cases.csv'), index=False)




if __name__ == "__main__":

    from simtools.Analysis.AnalyzeManager import AnalyzeManager
    from simtools.SetupParser import SetupParser

    from simulation.load_paths import load_box_paths
    data_path, project_path = load_box_paths(country_name='Burundi')

    SetupParser.default_block = 'HPC'
    SetupParser.init()

    working_dir = os.path.join(project_path, 'simulation_output', 'seasonality_calibration')
    start_year = 2015#2011
    end_year = 2018#2021

    expt_ids = {
        # 'TEST17blonger_2017_seasonality_calibration_1arch_Gitega_round2_iter0': '9aeeb00e-328d-eb11-a2ce-c4346bcb1550',
        'seasonality_calibration_1arch_Gitega_2017_round1_iter10': '979c5712-5b8d-eb11-a2ce-c4346bcb1550',

    }
    for expname, expid in expt_ids.items() :
        print('running expt %s' % expname)
        cur_monthlyTreatedCasesAnalyzer = monthlyTreatedCasesAnalyzer(expt_name=expname,
                                                                      channels=['Received_NMF_Treatment', 'Received_Treatment'],
                                                                      sweep_variables=["Run_Number", "__sample_index__"],
                                                                      working_dir=working_dir,
                                                                      start_year=start_year,
                                                                      end_year=end_year)

        analyzers = [cur_monthlyTreatedCasesAnalyzer]
        am = AnalyzeManager(expid, analyzers=analyzers, force_analyze=True)
        am.analyze()
