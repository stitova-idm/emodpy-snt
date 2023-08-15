import os
import manifest
import params
from snt.helpers_sim_setup import update_basic_params, set_input_files, set_drug_params


def set_config(config):
    update_basic_params(config, manifest, manifest.project_path)

    set_input_files(config, params.rep_admin, params.rep_admin,
                    population_size=params.simulation_pop)  # climate and demographics files

    set_drug_params(config)

    config.parameters.Simulation_Duration = (params.throwaway + 1) * 365

    if params.pull_from_serialization:
        from idmtools.core.context import get_current_platform
        platform = get_current_platform()

        ser_df = platform.create_sim_directory_df(params.burnin_id)
        ser_df['MaxHab'] = ser_df['MaxHab'].astype(float)
        if 'Admin_Name' in ser_df.columns.values:
            ser_df = ser_df[ser_df['MaxHab'] == params.burnin_max_habitat_value]
            ser_df = ser_df.set_index('Admin_Name')
            ser_path = ser_df.at[params.rep_admin, 'outpath']
        else:
            ser_path = ser_df['outpath'].values[0]

        config.parameters.Serialized_Population_Reading_Type = 'READ'
        config.parameters.Serialized_Population_Path = os.path.join(str(ser_path), 'output')
        config.parameters.Serialized_Population_Filenames = ['state-%05d.dtk' % (params.burnin_years * 365)]
        config.parameters.Enable_Random_Generator_From_Serialized_Population = 0
        config.parameters.Serialization_Mask_Node_Read = 16

    # move to here from set_input_files
    config.parameters.Demographics_Filenames = [params.demographics_file]

    return config
