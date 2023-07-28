# sometimes, the population size needs to change between calibration and later simulations (e.g., memory limitations)
#   the goal of this script is to update the demographics file and the vector monthly larval habitat file
#   to produce rescaled simulations with the desired human population size.

import pandas as pd
import numpy as np
import json
import os
import sys
sys.path.append('../../')
from simulation.load_paths import load_box_paths

data_path, project_path = load_box_paths(country_name='Example')


def update_demog_and_vector_pop(demog_filepath, larval_habitat_filepath, new_pop_size):
    # read in original files and get original population size
    orig_larval = pd.read_csv(larval_habitat_filepath)
    f = open(demog_filepath)
    demog_dict = json.load(f)
    original_pop_size = demog_dict['Nodes'][0]['NodeAttributes']['InitialPopulation']

    # calculate and update demographics parameters to correspond with new population size
    demog_death_rate = demog_dict['Defaults']['IndividualAttributes']['MortalityDistribution']['ResultValues'][0][0]
    demog_death_rate_scale = demog_dict['Defaults']['IndividualAttributes']['MortalityDistribution']['ResultScaleFactor']
    new_birth_number = new_pop_size * demog_death_rate_scale * demog_death_rate
    demog_dict_new = demog_dict
    demog_dict_new['Nodes'][0]['NodeAttributes']['InitialPopulation'] = new_pop_size
    demog_dict_new['Nodes'][0]['NodeAttributes']['BirthRate'] = new_birth_number

    # calculate and update larval habitat parameters to scale with new human population size
    new_larval = orig_larval
    add_larval_pop_factor = np.log10(new_pop_size / original_pop_size)
    new_larval['Constant'] = orig_larval['Constant'] + add_larval_pop_factor
    new_larval['MaxHab'] = orig_larval['MaxHab'] + add_larval_pop_factor

    # set new filenames and save files
    new_larval_filepath = larval_habitat_filepath.replace('.csv', '_%i.csv' % new_pop_size)
    new_larval.to_csv(new_larval_filepath, index=False)
    new_demog_filepath = demog_filepath.replace('.json', '_%i.json' % new_pop_size)
    # Serializing json and writing to .json file
    json_object = json.dumps(demog_dict_new, indent=4)
    with open(new_demog_filepath, "w") as outfile:
        outfile.write(json_object)

if __name__ == "__main__":
    # original demographics filepath
    demog_filepath = os.path.join(project_path, 'simulation_inputs', 'demographics_and_climate', '_entire_country', 'demographics_each_admin_10000.json')
    # original vector larval habitat filepath that was created with the original demographics file
    larval_habitat_filepath = os.path.join(project_path, 'simulation_inputs', 'larval_habitats', 'monthly_habitats_1.csv')
    # target demographic population size that will be used for later simulations
    new_pop_size = 6000
    # create new input files, rescaled for the new population size
    update_demog_and_vector_pop(demog_filepath, larval_habitat_filepath, new_pop_size)
