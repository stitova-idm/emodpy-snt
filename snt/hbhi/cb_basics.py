from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder

def initialize_cb (years, serialize, defaults='MALARIA_SIM', yr_plusone=True):
    cb = DTKConfigBuilder.from_defaults(defaults)
    cb.update_params({'Simulation_Duration': years*365+yr_plusone})
    if serialize :
        cb.update_params({
            'Serialization_Time_Steps' : [365*years],
            'Serialization_Type': 'TIMESTEP',
            'Serialization_Mask_Node_Write': 0,  
            # 0 corresponds to the previous version default: the same larval habitat 
            # parameters will be used in the burnin and pickup (from the burnin config)
            'Serialization_Precision': 'REDUCED'
        })
    else:
        cb.update_params({
            'Serialization_Type': 'NONE'
        })
