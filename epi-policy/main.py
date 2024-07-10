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
        x0 = [initial_scale_factor],

    )

    result = minimize(
        logistic_growth_objective,
        x0=[initial_scale_factor],  # initial guess as a list since L-BFGS-B is used for multidimensional optimization
        args=(y_max, mid_point, x_transition_units, x_vector),
        method='L-BFGS-B',
        bounds=bounds
    )



        