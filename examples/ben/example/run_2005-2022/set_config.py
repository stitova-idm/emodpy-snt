import manifest
import params
# from snt.helpers_sim_setup import update_basic_params
from emodpy_malaria.vector_config import set_species_param

from snt.hbhi.set_up_general import initialize_config


def set_config(config):
    # update_basic_params(config, params.project_path)

    # BASIC SETUP
    initialize_config(config, manifest, params.years, serialize=params.serialize)
    # cb.update_params({
    #     'x_temporary_Larval_Habitat': 1,  # Package default is 0.2
    #     'x_Base_Population': 1,
    #     'x_Birth': 1
    # })
    config.parameters.x_Temporary_Larval_Habitat = 1
    config.parameters.x_Base_Population = 1
    config.parameters.x_Birth = 1

    # cb.update_params({
    #     "Report_Event_Recorder": 1,
    #     "Report_Event_Recorder_Individual_Properties": [],
    #     "Report_Event_Recorder_Events": ['Received_NMF_Treatment', 'Received_Severe_Treatment',
    #                                      'Received_Treatment'],
    #     "Report_Event_Recorder_Ignore_Events_In_List": 0
    # })
    config.parameters.Report_Event_Recorder = 1
    config.parameters.Report_Event_Recorder_Individual_Properties = []
    config.parameters.Report_Event_Recorder_Events = ['Received_NMF_Treatment', 'Received_Severe_Treatment',
                                                      'Received_Treatment']
    config.parameters.Report_Event_Recorder_Ignore_Events_In_List = 0

    # update_species_param(cb, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
    # update_species_param(cb, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
    # update_species_param(cb, 'funestus', 'Anthropophily', 0.5, overwrite=True)
    # update_species_param(cb, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
    # update_species_param(cb, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
    # update_species_param(cb, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)

    set_species_param(config, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
    set_species_param(config, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Anthropophily', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
    set_species_param(config, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
    set_species_param(config, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)

    return config
