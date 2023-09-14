# emodpy-snt
Subnational(SNT) Workflow for 02_seasonality_calibration.py

We break one script into several smaller pieces and each piece focus only on one task.

manifest.py
- this script is used to set up locations for some important objects, such as Eradication, schema.json, etc.
- this file is required by emodpy-malaria

params.py
- this script is used to hold all defined parameters, 
- for example, user can just modify this file for different simulation sweeping or intervention sweeping

set_config.py
- This script is the place where user cab define/modify any configuration parameters

config_task.py
- This script is the where EMODTask got created (EMODTask is similar to DTKConfiguration in DTK-TOOLS)
- Here we connect all pieces (other smaller scripts) together
- Here user can add more information to task, such as adding reports...

run_calibration.py
- This is the entry-point to run calibration.

Note: only entry-point script run_calibration.py and manifest.py are "Must Have" in workflow, all other
scripts are optional (for a simple script, there is no need to break script into many pieces, we can
actually put all the coding into one script run_simulations.py)



