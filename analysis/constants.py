"""Global constants used in the prefactor fitting and for the analysis."""

import numpy as np

# Kappa from 10.1073/pnas.1917195117
kappa_0_spce = 5e11  # ns^-1 nm^-3
kappa_0_3D = kappa_0_spce
gamma_spce = 55.0  # spce water air surface tension in mN/m


kappa_0_2D_fit_params = np.array(
    [1.92962379e-35, 2.14682509e-01, 3.70256997e00, 1.26138372e01, 2.97690333e01]
)
kapp_0_2D_fit_sign = np.array([1, 1, -1, 1, -1])


def kappa_0_2D(contact_angle_deg):
    """Hetergenous cavitation prefactor.

    Values are hardcoded for the system of the puplication.

    Parameters
    ----------
    contact_angle_deg : float
        the contact angle in degrees"""
    kernel = 0
    for index, param in enumerate(kappa_0_2D_fit_params[1:]):
        exponent = index + 1
        kernel += kapp_0_2D_fit_sign[exponent] * (contact_angle_deg / param) ** exponent

    return kappa_0_2D_fit_params[0] * np.exp(kernel)
