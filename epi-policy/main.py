import numpy as np
import pandas as pd
from model import EpiModel

model = EpiModel("/Users/abhay/Documents/XLab/epi-policy/data/metapopulation-inputs-master.xlsx")

days = 365
model.run_simulation(days)