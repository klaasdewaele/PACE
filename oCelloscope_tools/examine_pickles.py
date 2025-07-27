#!/home/kdewaele/.conda/envs/ocelloscope_env/bin/python

import pandas as pd

path = "/home/kdewaele/ocelloscope/output/summary/performance_A_A_0.9_time_0-172800_bias_[2]_CA_0.9_2025_amfo_48u_18-07-2025-11-35-18.pkl"

df = pd.read_pickle(path)

print(df)