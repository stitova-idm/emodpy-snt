import os

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", "..", "..", ".."))
BEN_DIR = os.path.abspath(os.path.join(CURRENT_DIR, "..", ".."))

# change to your input folder which contains the required files...
input_dir = os.path.join(ROOT_DIR, "inputs")
download_dir = os.path.join(ROOT_DIR, "download")
schema_file = os.path.join(download_dir, "schema.json")
eradication_path = os.path.join(download_dir, "Eradication")

# SIF PATH
sif_path = '/projects/b1139/images/centos_dtk-build.sif'
