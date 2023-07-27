import emodpy_malaria.malaria_config as malaria_config
import manifest


def initialize_config(config, manifest, years, serialize, yr_plusone=True):
    config = malaria_config.set_team_defaults(config, manifest)
    config.parameters.Simulation_Duration = years * 365 + yr_plusone
    if serialize:
        if serialize:
            config.parameters.Serialization_Time_Steps = [365 * years]
            config.parameters.Serialization_Type = "TIMESTEP"
            config.parameters.Serialized_Population_Writing_Type = "TIMESTEP"
            config.parameters.Serialization_Mask_Node_Write = 0
            config.parameters.Serialization_Precision = "REDUCED"
    # else: nothing to do, we're already not serializing
