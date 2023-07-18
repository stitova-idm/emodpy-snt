import manifest
import params
from snt.helpers_sim_setup import update_basic_params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

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

    # move to here from set_input_files
    config.parameters.Demographics_Filenames = [params.demographics_file]

    return config
