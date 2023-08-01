import manifest
import params
from snt.helpers_sim_setup import update_basic_params
from emodpy_malaria import malaria_config

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

    # set DrugParams
    # Amodaquine
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Cmax", 270)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Decay_T1", 0.7)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Decay_T2", 15.9)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_PKPD_C50", 55)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Drug_Vd", 1)
    malaria_config.set_drug_param(config, 'Amodiaquine', "Max_Drug_IRBC_Kill", 0.2)

    # SulfadoxinePyrimethamine
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_Decay_T1", 11.5)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_Decay_T2", 11.5)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Drug_PKPD_C50", 0.9)
    malaria_config.set_drug_param(config, 'SulfadoxinePyrimethamine', "Max_Drug_IRBC_Kill", 0.28)

    return config
