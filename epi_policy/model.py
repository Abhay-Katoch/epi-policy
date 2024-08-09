import numpy as np
import pandas as pd
from epi_policy.functions import calibrate_logistic_function, logistic_growth_function

class EpiModel:

    def __init__(self, input_spreadsheet, generate_beta = True):

      self.inputs = pd.ExcelFile(input_spreadsheet)

      self.jurisdictions = self.inputs.parse("jurisdiction").set_index('jurisdiction.id')
      self.commuting_areas = self.inputs.parse("commuting_area")
      self.healthcosts = self.inputs.parse("healthcosts").set_index('severity')

      parameters = self.inputs.parse("parameters").set_index('parameter')
      parameter_dict = parameters["baseline"].to_dict()

      for key, value in parameter_dict.items():
        setattr(self, key, value)

      self.number_jurisdictions = len(self.jurisdictions.index)
      
      self.coordination = np.identity(len(self.jurisdictions.index))

    def generate_beta(self, work_travel_mixing = None):
      mixing_modes = ["Home", "Work", "Travel"]
      number_modes = len(mixing_modes)

      k_mat = np.array([self.k_home, self.k_work_travel, self.k_other] * self.number_jurisdictions)
      k_mat = k_mat.reshape(self.number_jurisdictions, number_modes).astype(float)

      # TODO: Implement try-catch routine to see if rows are not normalized

      home_mixing = np.diag(np.ones(self.number_jurisdictions))
      other_mixing = np.diag(np.ones(self.number_jurisdictions))

      if work_travel_mixing is None:
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

      self.beta = np.transpose(cbeta_i * np.transpose(mixing_matrix))

    def initialize_state(self, days):
        self.S = self.E = self.P = self.I = self.A = self.R = self.D = np.zeros(self.number_jurisdictions, dtype = int)
        self.N = np.array(self.jurisdictions["population"].values, dtype = int)

        self.survey_lag = np.full((self.number_jurisdictions), self.total_surv_lag, dtype = int)
        self.case_ascertainment = np.full((self.number_jurisdictions), self.p, dtype = float)
        self.incidence_history = np.zeros((days, self.number_jurisdictions), dtype = int)

        self.L_star_t = np.zeros((self.number_jurisdictions, 1))
        self.L_t = np.zeros((self.number_jurisdictions, 1))

        r_scale_factor = calibrate_logistic_function(
          y_max = self.p_r_star, 
          mid_point = self.t_mid_IFR, 
          x_transition_units = self.t_r_star, 
          x_vector = np.arange(0, 365)
        )

        IFR_time_mult = logistic_growth_function(
            y_max = self.p_r_star, 
            mid_point = self.t_mid_IFR, 
            x_transition_units = self.t_r_star, 
            x = np.arange(0, 365, 0.01), 
            scale_factor = r_scale_factor
        )

        self.time_varying_IFR = self.r * (1 - IFR_time_mult)
    
    def iterate_model(self, day):
        daily_incidences = self.incidence_history[(day - self.survey_lag), np.arange(self.number_jurisdictions)]
        C_hat_t = np.random.binomial(daily_incidences, self.case_ascertainment)
        x_star_t = np.minimum((1e5 * C_hat_t) / (self.N * self.case_ascertainment * self.c), self.L_max)

        self.L_star_t += 0.5 * (x_star_t[:, np.newaxis] - self.L_star_t)

        update_up = (day % self.a_up == 0) & (np.round(self.L_star_t) > self.L_t)
        self.L_t = np.where(update_up, np.floor(self.L_star_t), self.L_t)

        update_down = (day % self.a_down == 0) & (np.round(self.L_star_t) < self.L_t)
        self.L_t = np.where(update_down, np.floor(self.L_star_t), self.L_t)

        self.beta_t = (1 - (self.L_t * self.tau)) * self.beta
        self.lambda_t = np.matmul(self.beta_t, (self.P + self.I + self.A) / self.N)

        S_E = np.random.binomial(self.S, 1 - np.exp(-self.lambda_t))
        self.incidence_history[day] = S_E

        E_P = np.random.binomial(self.E, 1 - np.exp(-self.sigma))
        
        P_IA = np.random.binomial(self.P, 1 - np.exp(-self.delta))
        P_I = np.random.binomial(P_IA, 1 - self.rho)
        P_A = P_IA - P_I

        I_RD = np.random.binomial(self.I, 1 - np.exp(-self.gamma))
        I_D = np.random.binomial(I_RD, self.time_varying_IFR[day])
        I_R = I_RD - I_D

        A_R = np.random.binomial(self.A, 1 - np.exp(-self.gamma))

        ## Difference equations
        self.S = self.S - S_E
        self.E = self.E + S_E - E_P
        self.P = self.P + E_P - P_IA
        self.I = self.I + P_I - I_R
        self.A = self.A + P_A - A_R
        self.R = self.R + A_R + I_R
        self.D = self.D + I_D

    def run_simulation(self, days):
        self.S = self.jurisdictions["S0"].values.astype(int)
        self.I = self.jurisdictions["I0"].values.astype(int)

        columns = ['day', 'jurisdiction', 'NPI', 'S', 'E', 'P', 'I', 'A', 'R', 'D']
        self.results = pd.DataFrame(columns = columns)

        for day in range(days - 1):
          self.iterate_model(day)

          outputs = np.stack([
            np.repeat(day, self.number_jurisdictions), 
            range(self.number_jurisdictions), 
            self.L_t.flatten(), 
            self.S, self.E, self.P, self.I, self.A, self.R, self.D], axis=-1
          ).reshape(-1, 10)

          self.results = pd.concat([self.results, pd.DataFrame(outputs, columns = columns)])
      
