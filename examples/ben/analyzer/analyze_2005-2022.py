import argparse
import os
from simtools.Analysis.SSMTAnalysis import SSMTAnalysis
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from analyzer_collection import MonthlyPfPRAnalyzerU5, MonthlyTreatedCasesAnalyzer, AnnualAgebinPfPRAnalyzer,MonthlyAgebinPfPRAnalyzer, annualSevereTreatedByAgeAnalyzer
from load_paths import load_paths

SetupParser.default_block = 'NUCLUSTER'
iopath = load_paths()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-name', dest='exp_name', type=str, required=True)
    parser.add_argument('-id', dest='exp_id', type=str, required=True)

    return parser.parse_args()


working_dir = "."

if __name__ == "__main__":
    SetupParser.init()
    args = parse_args()

    experiments = {args.exp_name: args.exp_id}
    start_year = 2005
    end_year = 2023
    sweep_variables = ["Sample_ID",
                       "DS_Name",
                       'archetype',
                       'Run_Number']
    output_dir = os.path.join(iopath, 'simulation_output')

    for expt_name, expt_id in experiments.items():
        analyzer = [
            MonthlyPfPRAnalyzerU5(expt_name=expt_name,
                                  sweep_variables=sweep_variables,
                                  start_year=start_year,
                                  end_year=end_year,
                                  working_dir=os.path.join(output_dir, expt_name)),
            MonthlyTreatedCasesAnalyzer(expt_name=expt_name,
                                        sweep_variables=sweep_variables,
                                        working_dir=os.path.join(output_dir, expt_name),
                                        start_year=start_year,
                                        end_year=end_year),
            AnnualAgebinPfPRAnalyzer(expt_name=expt_name,
                                     sweep_variables=sweep_variables,
                                     working_dir=os.path.join(output_dir, expt_name),
                                     start_year=start_year,
                                     end_year=end_year-1),
            MonthlyAgebinPfPRAnalyzer(expt_name=expt_name,
                                     sweep_variables=sweep_variables,
                                     working_dir=os.path.join(output_dir, expt_name),
                                     start_year=start_year,
                                     end_year=end_year),
            annualSevereTreatedByAgeAnalyzer(expt_name=expt_name,
                                             sweep_variables=sweep_variables,
                                             working_dir=os.path.join(output_dir, expt_name),
                                             start_year=start_year,
                                             end_year=end_year-1)
        ]

        am = AnalyzeManager(expt_id, analyzers=analyzer)
        am.analyze()

