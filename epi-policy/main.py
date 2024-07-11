import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize

def logistic_growth_function(y_max, mid_point, x_transition_units, scale_factor, x):
    """
        Produces a logistic growth function given a specific value, x.

        Parameters
        ----------
        y_max : float
            The asymptote of the function, the maximum value that the function can reach.
        mid_point : float
            The x-value at which the function's output reaches half of the y_max.
        x_transition_units : float
            The scale of the transition, representing the range over which the non-asymptotic growth occurs.
        x : float or array_like
            The x-values at which the function needs to be evaluated.
        scale_factor : float, optional
            A scaling factor that adjusts the sharpness of the transition.

        Returns
        -------
        float or ndarray
            The output of the logistic function evaluated at x. Returns a single float if x is a scalar, or a numpy array if x is array-like.

        Examples
        --------
        >>> logistic_fn(100, 50, 10, np.array([45, 50, 55]))
        array([27.5 , 50.  , 72.5])
    """

    difference = mid_point - x
    scale_parameter = x_transition_units * scale_factor

    function = y_max / np.exp(difference * scale_parameter)

    return function

def logistic_growth_objective(scale_factor, y_max, mid_point, x_transition_units, x_vector):
    """
        Objective function for the calibration of the logistic growth function. This function 
        calculates the logistic growth values for given parameters and evaluates how well these 
        values match a desired transition, quantifying the difference through a squared error.

        Parameters
        ----------
        scale_factor: float
            Scale factor to be used in the logistic growth function.
        y_max : float
            Asymptote of the logistic growth function, representing the maximum achievable value.
        x_mid_point : float
            The x-value at which the logistic function reaches half of its maximum value.
        x_transition_units : float
            The scale of the transition, representing the range over which the non-asymptotic growth occurs.
        x_vector : array_like
            An array of x-values at which the logistic function needs to be evaluated.

        Returns
        -------
        float
            The calculated error between the observed transition range and the expected 
            transition duration. This value is used to assess the fit of the logistic model.
    """

    y_values = logistic_growth_function(y_max, mid_point, x_transition_units, scale_factor, x = x_vector)
    trans_range = np.where(np.diff(y_values) > 1e-6)[0]

    if len(trans_range) > 0:
        width = x_vector[trans_range[-1]] - x_vector[trans_range[0]]
    else:
        width = 0

    error = (x_transition_units - width) ** 2
    return error

def calibrate_logistic_function(y_max, mid_point, x_transition_units, x_vector):
    """
        Calibrates the logistic growth function by optimizing the scale factor to fit a specified transition range.

        Parameters
        ----------
        y_max : float
            Asymptote of the logistic growth function, the maximum value the function approaches as x increases.
        x_mid_point : float
            The x value at which the function value is half of y_max.
        x_trans : float
            Transition duration in terms of x values over which most of the growth occurs.
        x_vector : array_like, optional
            An array of x values over which the logistic function will be evaluated. If not provided,
            it defaults to a sequence from 0 to three times the x_mid_point, with a step of 0.01.

        Returns
        -------
        float
            The optimized scale factor that minimizes the discrepancy between the expected and
            calculated logistic function outputs over the specified x values.
    """

    initial_scale_factor = 4 / (x_transition_units ** 2)
    bounds = [(initial_scale_factor * 0.1, initial_scale_factor * 3)]

    result = minimize(
        logistic_growth_objective,
        x0=[initial_scale_factor],  # initial guess as a list since L-BFGS-B is used for multidimensional optimization
        args=(y_max, mid_point, x_transition_units, x_vector),
        method='L-BFGS-B',
        bounds=bounds
    )

    return result.x

class EpiModel:

    def __init__(self, input_spreadsheet): #**parameters):

        # Unpack dictionary of parameters into seperate variables of the form "self.key = value".
        #for key, value in parameters.items():
        #    setattr(self, key, value)

        # Alternative to the former if I refactor the code to simply access values from the dictionary of parameters.
        # self.parameters = parameters

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
      population_proportion = self.jurisdictions['population'].values / sum(self.jurisdictions['population'].values)

      c_1 = cbeta / (population_proportion[0] + sum(c_rr[1: ] * population_proportion[1: ]))
      cbeta_i = c_1 * c_rr

      k_i = cbeta_i * population_proportion / cbeta

      beta = np.transpose(cbeta_i * np.transpose(mixing_matrix))

      return beta

    def initialize_state(self):#, compartment_variables: list):

        #for variable in compartment_variables:
        #    setattr(self, variable, np.zeros((self.number_jurisdictions), dtype = int))

        #self.incidence_history = np.zeros((self.number_jurisdiction, self.total_survey_lag))
        #self.S = self.S0
        #self.I = self.I0

        self.incidence_history = np.zeros((self.number_jurisdictions, 15))
        self.lagged_incidence = np.zeros((self.number_jurisdictions))
        self.L = np.zeros((self.number_jurisdictions))
        self.NPI = np.zeros((self.number_jurisdictions))

        self.S = self.jurisdictions["S0"].values.astype(int)
        self.I = self.jurisdictions["I0"].values.astype(int)
        self.E = np.zeros((self.number_jurisdictions)).astype(int)
        self.P = np.zeros((self.number_jurisdictions)).astype(int)
        self.A = np.zeros((self.number_jurisdictions)).astype(int)
        self.R = np.zeros((self.number_jurisdictions)).astype(int)
        self.S_E = np.zeros((self.number_jurisdictions)).astype(int)
        self.E_P = np.zeros((self.number_jurisdictions)).astype(int)
        self.P_I = np.zeros((self.number_jurisdictions)).astype(int)
        self.P_IA = np.zeros((self.number_jurisdictions)).astype(int)
        self.P_A = np.zeros((self.number_jurisdictions)).astype(int)
        self.A_R = np.zeros((self.number_jurisdictions)).astype(int)
        self.I_R = np.zeros((self.number_jurisdictions)).astype(int)
        self.p_SE = np.zeros((self.number_jurisdictions)).astype(int)

        self.total_surv_lag = int(self.total_surv_lag)
        self.N = self.jurisdictions["population"].values

        column_names = ["day", "jurisdiction.id", "L", "NPI", "S", "E", "P", "I", "A", "R"]
        self.results = pd.DataFrame(columns = column_names)
        self.total_survey_lag = 13
        self.incidence_history = np.zeros((self.number_jurisdictions, self.total_survey_lag))

    def iterate_model(self, y, day):
        self.S, self.E, self.P, self.I, self.A, self.R = y.reshape(6,6)

        self.incidence_history[:, 1:] = self.incidence_history[:, :-1]
        self.incidence_history[:, 0] = np.random.binomial(self.S_E, self.p)

        lagged_incidence = self.incidence_history[:, self.total_survey_lag - 1]

        L_star_ind = np.minimum(1e5 * (lagged_incidence / self.N) / (self.c * self.p), self.L_max + 0.01) # The 0.01 is there just to floor the value
        L_star_matrix = np.maximum.accumulate(np.transpose(L_star_ind) * self.coordination, axis = 1)
        L_star_f = np.where(self.L_c, L_star_matrix[:, -1], L_star_ind)

        self.L = (L_star_f - self.L) / 2
        change_NPI = np.where(self.L > self.NPI, day % self.a_up == 0, day % self.a_down == 0)
        self.NPI = np.where(change_NPI, np.floor(self.L), self.NPI)

        variant_rr = 1# self.variant_beta_rr if day >= self.variant_t else 1
        total_infectious = self.P + self.I + self.A
        transmission_rate = self.beta * ((1 - self.NPI[:, np.newaxis] * self.tau) * variant_rr)
        lambda_vector = np.sum(transmission_rate * (total_infectious / self.N), axis = 1)

        self.S_E = np.random.binomial(self.S.astype(int), 1 - np.exp(-lambda_vector))
        self.E_P = np.random.binomial(self.E.astype(int), 1 - np.exp(-self.sigma))

        self.P_IA = np.random.binomial(self.P.astype(int), 1 - np.exp(-self.delta))
        self.P_I = np.random.binomial(self.P_IA.astype(int), 1 - self.rho)
        self.P_A = self.P_IA - self.P_I

        # Logic error in original implementation of I/A -> R; D is sampled from R, therefore there is a chance that D is sampled from A.
        #   - Unless r is adjusted for *all* cases and not just I, this leads to an oversampling of deaths.
        #   - Also, gamma should probably be different between I and A; it takes less time (assumed) to recover from A than I.
        self.I_R = np.random.binomial(self.I.astype(int), 1 - np.exp(-self.gamma))
        self.A_R = np.random.binomial(self.A.astype(int), 1 - np.exp(-self.gamma))

        dS = self.S - self.S_E
        dE = self.E + self.S_E - self.E_P
        dP = self.P + self.E_P - self.P_IA
        dI = self.I + self.P_I - self.I_R
        dA = self.A + self.P_A - self.A_R
        dR = self.R + self.I_R + self.A_R

        return np.array([dS, dE, dP, dI, dA, dR]).flatten()

    def run_simulation(self, days):
        self.initialize_state()

        timescale = np.linspace(0, days, days)
        initial_values = np.array([self.S, self.E, self.P, self.I, self.A, self.R]).flatten().astype(int)
        solution = odeint(self.iterate_model, initial_values, timescale) # Pass the array of time points

        return solution