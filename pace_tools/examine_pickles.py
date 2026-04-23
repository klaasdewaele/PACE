#!/usr/bin/env python3
"""
examine_pickles.py — utility for inspecting pipeline result pickle files.

Reads a performance summary pickle file produced by set_thresh.py and prints
the contents as a DataFrame.

Usage:
    python examine_pickles.py <path_to_performance_pickle.pkl>
"""

import sys
import pandas as pd

if len(sys.argv) != 2:
    print(f"Usage: {sys.argv[0]} <path_to_performance_pickle.pkl>")
    sys.exit(1)

path = sys.argv[1]
df = pd.read_pickle(path)
print(df)
