import os
import pandas as pd
from snt.hbhi.set_up_general import load_master_csv, load_rel_abund_df
from emodpy_malaria.interventions.drug_campaign import add_drug_campaign
from emodpy_malaria.reporters.builtin import add_malaria_summary_report


def add_nmf_trt(campaign, num_years, start_year,
                coverage=0.0038, diagnostic_threshold=5):
    for nmf_years in range(num_years):
        add_drug_campaign(campaign, 'MSAT', 'AL',
                          start_days=[1 + 365 * nmf_years + start_year * 365],
                          coverage=coverage,
                          repetitions=365, tsteps_btwn_repetitions=1,
                          diagnostic_type='PF_HRP2', diagnostic_threshold=diagnostic_threshold,
                          receiving_drugs_event_name='Received_NMF_Treatment')


def tryread_df(csv_path):
    try:
        df = pd.read_csv(csv_path)
    except IOError:
        if not csv_path == '':
            print(f"WARNING: Cannot read {csv_path}.")
        df = pd.DataFrame()
    return df


def generate_multiyr_df(df, num_year, simday_col='simday'):
    dfs = []
    for year in range(num_year):
        df_yr = df.copy()
        df_yr[simday_col] = df_yr[simday_col] + 365 * year
        dfs.append(df_yr)
    return pd.concat(dfs)


def add_monthly_parasitemia_rep_by_year(task, manifest, num_year, tot_year, sim_start_year,
                                        yr_plusone=True, age_bins=None,
                                        prefix='Monthly_', ipfilter=''):
    if not age_bins:
        age_bins = [0.25, 5, 15, 30, 50, 125]
    if (num_year > tot_year):
        raise ValueError('num_year must be <= tot_year')
    for year in range(num_year):
        sim_year = tot_year - num_year + year
        start = yr_plusone + 365 * sim_year
        desc = prefix + '%d' % (sim_year + sim_start_year)
        add_malaria_summary_report(task, manifest, start_day=start, age_bins=age_bins,
                                   reporting_interval=30,
                                   end_day=start + 365, filename_suffix=desc,
                                   parasitemia_bins=[10, 50, 1e9],
                                   must_have_ip_key_value=ipfilter)


def add_annual_parasitemia_rep(task, manifest, num_year, tot_year, sim_start_year,
                               yr_plusone=True, age_bins=None,
                               prefix='Annual_', ipfilter=''):
    if not age_bins:
        age_bins = [0.25, 5, 15, 30, 50, 125]
    if (num_year > tot_year):
        raise ValueError('num_year must be <= tot_year')
    sim_year = tot_year - num_year
    start = yr_plusone + 365 * sim_year
    desc = prefix + '%dto%d' % (sim_year + sim_start_year, tot_year + sim_start_year - 1)
    add_malaria_summary_report(task, manifest, start_day=start, age_bins=age_bins,
                               reporting_interval=365, filename_suffix=desc,
                               parasitemia_bins=[10, 50, 1e9],
                               must_have_ip_key_value=ipfilter)


def read_main_dfs(projectpath, master_file=None, country=None,
                  larval_hab_csv='simulation_inputs/monthly_habitats.csv'):
    master_df = load_master_csv(projectpath=projectpath,
                                file=master_file,
                                country=country)
    rel_abund_df = load_rel_abund_df(projectpath)
    lhdf = pd.read_csv(os.path.join(projectpath, larval_hab_csv))

    return master_df, rel_abund_df, lhdf
