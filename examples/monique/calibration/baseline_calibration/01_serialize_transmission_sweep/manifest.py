import os
import params

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "..", "..", ".."))

# change to your input folder which contains the required files...
download_dir = os.path.join(PROJECT_DIR, "download")
schema_file = os.path.join(download_dir, "schema.json")
eradication_path = os.path.join(download_dir, "Eradication")

sif_id = None
