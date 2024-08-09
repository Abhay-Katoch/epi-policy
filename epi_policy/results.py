import numpy as np
import pandas as pd
from model import EpiModel

""" # Instantiate the model.
model = EpiModel("/Users/abhay/Documents/XLab/epi-policy/data/metapopulation-inputs-master.xlsx")

# Create the matrix with jurisdiction information: county name, population, and S0 + I0
pop = pd.read_csv("/Users/abhay/Documents/XLab/epi-policy/data/co-est2023-alldata.csv", encoding="latin-1")
mass = pop.loc[(pop["STATE"] == 25) & (pop["COUNTY"] != 0)]
n = len(mass)
proportion_infected = 0.00001

model.jurisdictions = pd.DataFrame({
    'jurisdiction.name': mass["CTYNAME"].values,  # Empty string values represented as NaN
    'npi.coordination.block': [1] * n,
    'commuting.area.id': [1] * n,
    'population': mass["ESTIMATESBASE2020"].values.astype(int),
    'S0': np.round(mass["ESTIMATESBASE2020"].values.astype(int) * (1 - proportion_infected)),
    'I0': np.round(mass["ESTIMATESBASE2020"].values.astype(int) * proportion_infected),  # Pre-filled with values from your example
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
model.save_results(file_location = "/Users/abhay/Documents/XLab/epi-policy/results/raw_results_long.csv") """

# Figure 1 - Marginal Benefit of ESS w/o heterogeneity
np.random.seed(0)

days = 365
outcomes = ["deaths_per_100k", "CH_illness", "CH_deaths", "CH", "CNPI", "C"]
marginal_NMB_results = pd.DataFrame(columns = outcomes)

base_model = EpiModel("/Users/abhay/Documents/XLab/epi-policy/data/metapopulation-inputs-master.xlsx")
base_model.generate_beta(work_travel_mixing = None)
base_model.initialize_state(days)
base_model.run_simulation()
base_model.save_results()

marginal_NMB_results.loc["Base"] = base_model.results.iloc[-1][outcomes]

model_no_het = EpiModel("/Users/abhay/Documents/XLab/epi-policy/data/metapopulation-inputs-master.xlsx")
model_no_het.generate_beta(work_travel_mixing = None)
model_no_het.initialize_state(days)

model_no_het.survey_lag = 8
model_no_het.case_ascertainment = 1

model_no_het.run_simulation()
model_no_het.save_results()
marginal_NMB_results.loc[f"w/ ESS"] = model_no_het.results.iloc[-1][outcomes]


"""
for jurisdiction, name in enumerate(model_no_het.jurisdictions["jurisdiction.name"]):
    model_no_het.initialize_state(days)
    model_no_het.survey_lag[:jurisdiction + 1] = 3
    model_no_het.case_ascertainment[:jurisdiction + 1] = 1

    model_no_het.run_simulation()
    model_no_het.save_results()

    marginal_NMB_results.loc[f"{jurisdiction + 1} w/ ESS"] = model_no_het.results.iloc[-1][outcomes]
"""