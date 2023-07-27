import os
import pandas as pd
from simulation.load_paths import load_box_paths

data_path, project_path = load_box_paths(country_name='Burundi')

wdir = os.path.join(project_path, 'simulation_inputs')
future_projections = True

while True :
    if future_projections:
        scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_for_projections.csv')
    else:
        scenario_fname = os.path.join(project_path, 'simulation_inputs', '_intervention_file_references', 'Interventions_to_present.csv')
    df = pd.read_csv(scenario_fname)

    if len(df[df['status'] == 'run']) == 0:
        break

    if future_projections:
        os.system('python ./run_future_scenarios.py')
    else:
        os.system('python ./run_to_present.py')

    i = df[df['status'] == 'run'].index[0]
    df.loc[i, 'status'] = 'queued'
    df.to_csv(scenario_fname, index=False)

