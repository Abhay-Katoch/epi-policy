import numpy as np
import pandas as pd

class EpiModel:

    def __init__(self, input_spreadsheet):

      self.inputs = pd.ExcelFile(input_spreadsheet)

      self.jurisdictions = self.inputs.parse("jurisdiction").set_index('jurisdiction.id')
      self.commuting_areas = self.inputs.parse("commuting_area")
      self.healthcosts = self.inputs.parse("healthcosts").set_index('severity')

      parameters = self.inputs.parse("parameters").set_index('parameter')
      parameter_dict = parameters["baseline"].to_dict()

      for key, value in parameter_dict.items():
            setattr(self, key, value)

      self.number_jurisdictions = len(self.jurisdictions.index)

      self.beta = self.generate_beta()
      self.coordination = np.identity(len(self.jurisdictions.index))

    def generate_beta(self):
      mixing_modes = ["Home", "Work", "Travel"]
      number_modes = len(mixing_modes)

      k_mat = np.array([self.k_home, self.k_work_travel, self.k_other] * self.number_jurisdictions)
      k_mat = k_mat.reshape(self.number_jurisdictions, number_modes).astype(float)

      # TODO: Implement try-catch routine to see if rows are not normalized

      home_mixing = np.diag(np.ones(self.number_jurisdictions))
      other_mixing = np.diag(np.ones(self.number_jurisdictions))
      work_travel_mixing = np.zeros((self.number_jurisdictions, self.number_jurisdictions))

      for i in self.jurisdictions.index:
          id_i = self.jurisdictions.loc[self.jurisdictions.index == i, 'commuting.area.id'].values[0]

          intra_rate = self.commuting_areas.loc[self.commuting_areas['commuting.area.id'] == id_i, 'intra.commuting.rate'].values[0]
          inter_rate = self.commuting_areas.loc[self.commuting_areas['commuting.area.id'] == id_i, 'inter.commuting.rate'].values[0]

          for j in self.jurisdictions.index:
              # Check if j is in the same commuting area as i
              id_j = self.jurisdictions.loc[self.jurisdictions.index == j, 'commuting.area.id'].values[0]
              is_same_commuting_area = (id_i == id_j)

              if (i != j):
                if (is_same_commuting_area):
                  work_travel_mixing[i - 1, j - 1] = intra_rate
                else:
                  work_travel_mixing[i - 1, j - 1] = inter_rate

          work_travel_mixing[i - 1, i - 1] = 1 - sum(work_travel_mixing[i - 1, ])

      M = np.stack([home_mixing, work_travel_mixing, other_mixing])

      mixing_matrix = np.zeros((self.number_jurisdictions, self.number_jurisdictions))

      for i, mode in enumerate(mixing_modes):
          mixing_matrix = mixing_matrix + np.multiply(k_mat[:, i], M[i])

      c_rr = self.jurisdictions['mixing.risk.ratio'].astype(float).values

      tau_eff = (1 / self.delta) + (1 / self.gamma)
      cbeta = self.R0 / tau_eff
      population_proportion = self.jurisdictions['population'].values.astype(int) / sum(self.jurisdictions['population'].values)

      c_1 = cbeta / (population_proportion[0] + sum(c_rr[1: ] * population_proportion[1: ]))
      cbeta_i = c_1 * c_rr

      k_i = cbeta_i * population_proportion / cbeta

      beta = np.transpose(cbeta_i * np.transpose(mixing_matrix))

      return beta

    def initialize_state(self, days):#, compartment_variables: list):

        self.S = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.I = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.E = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.P = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.A = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.R = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.D = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.S_E = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.E_P = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.P_I = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.P_A = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.A_R = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.I_R = np.zeros((days, self.number_jurisdictions)).astype(int)
        self.I_D = np.zeros((days, self.number_jurisdictions)).astype(int)

        self.L_star_t = np.zeros((self.number_jurisdictions))
        self.L_t = np.zeros((self.number_jurisdictions))
        self.N = self.jurisdictions["population"].values
        self.total_survey_lag = 13
    
    def iterate_model(self, day):
        C_hat_t = np.random.binomial(self.S_E[day - self.total_survey_lag], self.p)
        x_star_t = np.minimum((1e3 * C_hat_t) / (self.N * self.p * self.c), self.L_max)

        self.L_star_t += 0.5 * (x_star_t - self.L_star_t)

        update_up = (day % self.a_up == 0) & (np.round(self.L_star_t) > self.L_t)
        self.L_t = np.where(update_up, np.round(self.L_star_t), self.L_t)

        update_down = (day % self.a_down == 0) & (np.round(self.L_star_t) < self.L_t)
        self.L_t = np.where(update_down, np.round(self.L_star_t), self.L_t)

        self.beta_t = np.dot((1 - (self.L_t * self.tau)), self.beta)
        self.lambda_t = self.beta_t * (self.P[day] + self.I[day] + self.A[day])

        self.S_E[day] = np.random.binomial(self.S[day].astype(int), 1 - np.exp(- self.lambda_t / self.N))
        self.E_P[day] = np.random.binomial(self.E[day].astype(int), 1 - np.exp(- self.sigma))
        
        # This may be an error; if I and A are independently sampled from P, then there is a possibility that I_t + A_t > P_{t - 1}
        self.P_I[day] = np.random.binomial(self.P[day].astype(int), 1 - np.exp(- self.delta * (1 - self.rho)))
        self.P_A[day] = np.random.binomial(self.P[day].astype(int), 1 - np.exp(- self.delta * self.rho))

        # Logic error in original implementation of I/A -> R; D is sampled from R, therefore there is a chance that D is sampled from A.
        #   - Unless r is adjusted for *all* cases and not just I, this leads to an oversampling of deaths.
        #   - Also, gamma should probably be different between I and A; it takes less time (assumed) to recover from A than I.
        self.I_R[day] = np.random.binomial(self.I[day].astype(int), 1 - np.exp(-self.gamma  * (1 - self.r)))
        self.A_R[day] = np.random.binomial(self.A[day].astype(int), 1 - np.exp(-self.gamma))

        self.I_D[day] = np.random.binomial(self.I[day].astype(int), 1 - np.exp(-self.gamma  * (self.r)))

        self.S[day + 1] = self.S[day] - self.S_E[day]
        self.E[day + 1] = self.E[day] + self.S_E[day] - self.E_P[day]
        self.P[day + 1] = self.P[day] + self.E_P[day] - (self.P_I[day] + self.P_A[day])
        self.I[day + 1] = self.I[day] + self.P_I[day] - self.I_R[day]
        self.A[day + 1] = self.A[day] + self.P_A[day] - self.A_R[day]
        self.R[day + 1] = self.R[day] + self.I_R[day] + self.A_R[day]
        self.D[day + 1] = self.D[day] + self.I_D[day]

    def run_simulation(self, days: int):
        self.initialize_state(days)
        self.S[0] = self.jurisdictions["S0"].values.astype(int)
        self.I[0] = self.jurisdictions["I0"].values.astype(int)

        for day in range(days - 1):
            self.iterate_model(day)