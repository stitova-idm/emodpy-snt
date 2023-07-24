import os
import numpy as np
import pandas as pd
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer


class MonthlyPfPRU5Analyzer(BaseAnalyzer):

    def __init__(self, expt_name, sweep_variables=None, working_dir=".", start_year=2010, end_year=2016,
                 input_filename_base='MalariaSummaryReport_Monthly'):
        super(MonthlyPfPRU5Analyzer, self).__init__(working_dir=working_dir,
                                                    filenames=["output/%s%d.json" % (input_filename_base, x)
                                                               for x in range(start_year, (1+end_year))]
                                                    )
        self.sweep_variables = sweep_variables or ["Run_Number"]
        self.mult_param = 'Habitat_Multiplier'
        self.expt_name = expt_name
        self.start_year = start_year
        self.end_year = end_year


    def select_simulation_data(self, data, simulation):

        adf = pd.DataFrame()
        aa = 1  # index for the 0.5-5 year age group
        for year, fname in zip(range(self.start_year, (self.end_year+1)), self.filenames):
            # population size
            pop = data[fname]['DataByTimeAndAgeBins']['Average Population by Age Bin'][:12]  # remove final five days: assume final five days have same average as rest of month
            pop_monthly = [x[aa] for x in pop]
            # PfPR
            d = data[fname]['DataByTimeAndAgeBins']['PfPR by Age Bin'][:12]  # remove final five days: assume final five days have same average as rest of month
            pfpr_monthly = [x[aa] for x in d]
            # combine in data frame
            simdata = pd.DataFrame({'month': list(range(1, 13)),
                                    'Pop': pop_monthly,
                                    'PfPR U5': pfpr_monthly,
                                    })
            simdata['year'] = year
            adf = pd.concat([adf, simdata])

        for sweep_var in self.sweep_variables:
            if sweep_var in simulation.tags.keys():
                adf[sweep_var] = simulation.tags[sweep_var]
        return adf

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        if not os.path.exists(self.working_dir):
            os.mkdir(self.working_dir)

        all_df = pd.concat(selected).reset_index(drop=True)

        # don't take average across runs; calculate the best xLH for each seed (each of which uses different quantile)
        # all_df = all_df.groupby([self.mult_param, 'month', 'admin_name'])[['PfPR U5', 'Pop']].agg(np.mean).reset_index()
        # all_df = all_df.sort_values(by=[self.mult_param, 'month', 'archetype', 'DS_Name_for_ITN'])
        all_df.to_csv(os.path.join(self.working_dir, 'monthly_U5_PfPR.csv'), index=False)
