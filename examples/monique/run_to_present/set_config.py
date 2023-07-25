import manifest
import params
from snt.helpers_sim_setup import update_basic_params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

    config.parameters.Report_Event_Recorder = 1
    config.parameters.Report_Event_Recorder_Events = ['Received_Severe_Treatment']  # !, 'Bednet_Got_New_One'],
    config.parameters.Report_Event_Recorder_Ignore_Events_In_List = 0
    config.parameters.Report_Event_Recorder_Individual_Properties = []

    if params.serialize:
        config.parameters.Simulation_Duration = params.years * 365 + 1
        config.parameters.Serialization_Time_Steps = [365 * (params.years - 5), 365 * params.years]
        config.parameters.Serialized_Population_Reading_Type = 'NONE'
        config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        config.parameters.Serialization_Mask_Node_Write = 0
        config.parameters.Serialization_Precision = 'REDUCED'
    else:
        config.parameters.Serialized_Population_Reading_Type = 'NONE'
        config.parameters.Simulation_Duration = params.years * 365 + 1
        config.parameters.Serialized_Population_Writing_Type = 'NONE'

    return config
