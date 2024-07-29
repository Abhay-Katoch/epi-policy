import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from model import EpiModel

# Instantiate the model.
model = EpiModel("/Users/abhay/Documents/XLab/epi-policy/data/metapopulation-inputs-master.xlsx")

# Create the matrix with jurisdiction information: county name, population, and S0 + I0
pop = pd.read_csv("/Users/abhay/Documents/XLab/epi-policy/data/co-est2023-alldata.csv", encoding="latin-1")
mass = pop.loc[(pop["STATE"] == 25) & (pop["COUNTY"] != 0)]
n = len(mass)

model.jurisdictions = pd.DataFrame({
    'jurisdiction.name': mass["CTYNAME"].values,  # Empty string values represented as NaN
    'npi.coordination.block': [1] * n,
    'commuting.area.id': [1] * n,
    'population': mass["ESTIMATESBASE2020"].values.astype(int),
    'S0': mass["ESTIMATESBASE2020"].values.astype(int) - 200,
    'I0': [200] * n,  # Pre-filled with values from your example
    'cost.npi': [0.0] * n,  # Pre-filled with values from your example
    'mixing.risk.ratio': [1.0] * n  # Pre-filled with values from your example
}, index=pd.Index(range(n), name='jurisdiction.id'))

model.number_jurisdictions = len(model.jurisdictions)

# Make the commuting matrix from the proportion of commuters <-> each jurisdiction
ACS_flows = pd.read_excel("/Users/abhay/Documents/XLab/epi-policy/data/2016-2020 5-Year ACS Commuting Flows.xlsx", header = 7)
ACS_flows_mass = ACS_flows[(ACS_flows['State FIPS Code'] == "25") & (ACS_flows['State FIPS Code.1'] == 25)]

commuting_matrix = pd.pivot_table(ACS_flows_mass, 
                                  index='County Name', 
                                  columns='County Name.1', 
                                  values='Workers in Commuting Flow', 
                                  aggfunc='sum', 
                                  fill_value=0) 
total_commuters_per_origin = commuting_matrix.sum(axis=1)
beta = commuting_matrix.div(total_commuters_per_origin, axis='index')

model.generate_beta(work_travel_mixing=beta)

# Run Simulation
days = 365
model.run_simulation(days)

# Save model outputs to a .csv file. 
model_outputs = np.stack([model.S, model.E, model.P, model.I, model.A, model.R, model.D], axis=-1).reshape(-1, 7)

results_long = pd.DataFrame(model_outputs, columns=['S', 'E', 'P', 'I', 'A', 'R', 'D'])
results_long['day'] = np.repeat(np.arange(model.days), model.number_jurisdictions)
results_long['jurisdiction'] = np.tile(np.arange(model.number_jurisdictions), model.days)

results_long = results_long[['day', 'jurisdiction', 'S', 'E', 'P', 'I', 'A', 'R', 'D']]
results_long.to_csv("/Users/abhay/Documents/XLab/epi-policy/results/raw_results_long.csv")