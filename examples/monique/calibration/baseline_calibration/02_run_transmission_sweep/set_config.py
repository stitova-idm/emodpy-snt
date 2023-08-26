import manifest
import params
from snt.helpers_sim_setup import update_basic_params, set_drug_params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

    set_drug_params(config)

    config.parameters.Simulation_Duration = params.years * 365
    config.parameters.Serialized_Population_Writing_Type = 'NONE'  # TODO: remove

    # move to here from set_input_files
    config.parameters.Demographics_Filenames = [params.demographics_file]

    # add Custom Individual Events
    config.parameters.Custom_Individual_Events.extend(["Bednet_Got_New_One", "Bednet_Using", "Bednet_Discarded"])

    return config
