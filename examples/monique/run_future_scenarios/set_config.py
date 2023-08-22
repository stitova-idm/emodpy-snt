import manifest
import params
from snt.helpers_sim_setup import update_basic_params, set_drug_params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

    set_drug_params(config)

    config.parameters.Report_Event_Recorder = 1
    config.parameters.Report_Event_Recorder_Events = ['Received_Severe_Treatment', 'Bednet_Got_New_One']
    config.parameters.Report_Event_Recorder_Ignore_Events_In_List = 0
    config.parameters.Report_Event_Recorder_Individual_Properties = []

    if params.serialize:
        config.parameters.Simulation_Duration = params.years * 365 + 1
        config.parameters.Serialization_Time_Steps = [365 * (params.years - 5), 365 * params.years]
        config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        config.parameters.Serialization_Mask_Node_Write = 0
        config.parameters.Serialization_Precision = 'REDUCED'
    else:
        config.parameters.Simulation_Duration = params.years * 365 + 1

    # move to here from set_input_files
    config.parameters.Demographics_Filenames = [params.demographics_file]

    # add Custom Individual Events
    config.parameters.Custom_Individual_Events.extend(
        ['ReceivedTreatment', 'Received_Severe_Treatment', 'Received_NMF_Treatment'])

    return config
