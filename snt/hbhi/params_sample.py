import pandas as pd
from pyDOE import *
from scipy.stats.distributions import norm
from scipy.stats import uniform
from random import randint


def make_distribution(distro: str, param1: float, param2: float):
    if 'uniform' in distro:
        return uniform(loc=param1, scale=param2 - param1)
    elif 'normal' in distro:
        return norm(loc=param1, scale=param2)
    else:
        raise Exception("Sorry, not supporting other distribution yet.")


def transform(arr, method: str):
    if method == 'invlogit':
        return 1 / (1 + np.exp(-arr))
    elif method == 'exp':
        return np.exp(arr)
    elif method == 'exp10':
        return 10**arr
    else:
        raise Exception("Sorry not supporting other transformation yet.")


def gen_samples_from_df(df: pd.DataFrame, my_ds: str, nsamples: int = 1000):
    df1 = df[df.DS_Name == my_ds]
    design = lhs(len(df1), samples=nsamples)

    i = 0
    columns = []
    for r, row in df1.iterrows():
        distro = row['Distribution']
        param1 = row['Param1']
        param2 = row['Param2']

        if distro == 'normal' and param2 == 0:
            design[:, i] = param1
        else:
            design[:, i] = make_distribution(distro, param1, param2).ppf(design[:, i])
        design[:, i] = transform(design[:, i], row['Transform'])

        i += 1
        columns = columns + [row['Var_Name']]

    sample_df = pd.DataFrame(design, columns=columns)
    sample_df['DS_Name'] = my_ds
    sample_df['id'] = range(nsamples)
    sample_df['seed1'] = [randint(0, 60000) for j in range(nsamples)]
    sample_df['seed2'] = [randint(0, 60000) for j in range(nsamples)]
    return sample_df


def reduce_dim(df: pd.DataFrame, col: str, scale: int = 5, reduce_seed: bool = True):
    df1 = df.sort_values(col).reset_index(drop=True)
    df1.id = sorted(df1.id)
    df1['reduce_id'] = df1.id // scale
    df1 = df1.assign(tmp = lambda x: x.groupby('reduce_id').transform('max')[col])
    df1[col] = df1.tmp
    df1 = df1.drop(columns='tmp')
    if reduce_seed:
        df1 = df1.assign(seed1 = lambda x: x.groupby('reduce_id').transform('mean')['seed1'])
        df1.seed1 = df1.seed1.astype(int)

    return df1
