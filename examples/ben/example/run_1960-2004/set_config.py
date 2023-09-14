import manifest
import params
from emodpy_malaria.vector_config import set_species_param
from emodpy_malaria.malaria_config import set_drug_param
from snt.hbhi.set_up_general import initialize_config


def set_config(config):
    initialize_config(config, manifest, params.years, params.serialize)

    config.parameters.x_Temporary_Larval_Habitat = 1  # Package default is 0.2
    config.parameters.x_Base_Population = 1
    config.parameters.x_Birth = 1

    set_species_param(config, 'arabiensis', 'Anthropophily', 0.88, overwrite=True)
    set_species_param(config, 'arabiensis', 'Indoor_Feeding_Fraction', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Anthropophily', 0.5, overwrite=True)
    set_species_param(config, 'funestus', 'Indoor_Feeding_Fraction', 0.86, overwrite=True)
    set_species_param(config, 'gambiae', 'Anthropophily', 0.74, overwrite=True)
    set_species_param(config, 'gambiae', 'Indoor_Feeding_Fraction', 0.9, overwrite=True)

    # Update drug params
    set_drug_param(config, 'Amodiaquine', "Drug_PKPD_C50", 303)
    set_drug_param(config, 'Amodiaquine', "Max_Drug_IRBC_Kill", 3.76)
    set_drug_param(config, 'Amodiaquine', "Drug_Cmax", 95)
    set_drug_param(config, 'Amodiaquine', "Drug_Decay_T1", 0.775)
    set_drug_param(config, 'Amodiaquine', "Drug_Decay_T2", 37.5)
    set_drug_param(config, 'Amodiaquine', "Drug_Vd", 15.6)
    set_drug_param(config, 'Amodiaquine', "Fractional_Dose_By_Upper_Age",
                                  [{"Upper_Age_In_Years": 1, "Fraction_Of_Adult_Dose": 0.376},
                                   {"Upper_Age_In_Years": 5, "Fraction_Of_Adult_Dose": 0.44}])

    return config
