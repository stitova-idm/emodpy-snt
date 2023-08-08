import copy
import emod_api.config.default_from_schema_no_validation as dfs
from emodpy_malaria.malaria_config import configure_linear_spline, set_species_param
from snt.calibration.helpers_seasonality_calibration import get_burnin_spline_values
from snt.helpers_sim_setup import update_basic_params, set_input_files, set_drug_params

import manifest
import params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

    set_input_files(config, params.rep_admin, params.rep_admin,
                    population_size=params.simulation_pop)  # climate and demographics files

    set_drug_params(config)

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

    # set_input_files(cb, rep_admin, rep_admin, population_size=simulation_pop)
    #############################################################################################################
    # specify vectors and monthly habitat
    month_vals, fractions = get_burnin_spline_values(params.rep_admin, params.month_scalar, manifest.project_path)

    habitat = dfs.schema_to_config_subnode(manifest.schema_file, ["idmTypes", "idmType:VectorHabitat"])
    for (s, sp) in zip(fractions, ['arabiensis', 'funestus', 'gambiae']):
        linear_spline_habitat = configure_linear_spline(manifest,
                                                        max_larval_capacity=pow(10, params.max_habitat_value),
                                                        capacity_distribution_number_of_years=1,
                                                        capacity_distribution_over_time={
                                                            "Times": [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304,
                                                                      334],
                                                            "Values": month_vals
                                                        }
                                                        )

        set_species_param(config, sp, "Habitats", linear_spline_habitat, overwrite=True)

        new_habitat = copy.deepcopy(habitat)
        new_habitat.parameters.Habitat_Type = "CONSTANT"
        new_habitat.parameters.Max_Larval_Capacity = pow(10, params.max_habitat_value) * s
        set_species_param(config, sp, "Habitats", new_habitat.parameters)

    return config
