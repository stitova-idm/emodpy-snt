import os
import params

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "..", "..", ".."))

# change to your input folder which contains the required files...
input_dir = os.path.join(params.project_path, "simulation_inputs")
download_dir = os.path.join(PROJECT_DIR, "download")
schema_file = os.path.join(download_dir, "schema.json")
eradication_path = os.path.join(download_dir, "Eradication")

relative_path = os.path.join('demographics_and_climate', '_entire_country')
demog_file = f'demographics_each_admin_{params.population_size}.json'

plugins_folder = os.path.join(download_dir, "plugins")


sif_id = None
