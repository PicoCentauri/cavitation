import numpy as np
import scipy
from constants import kappa_0_2D, kappa_0_3D, kappa_0_spce
from prefactor_fit import cavitation_rate


def waiting_time(contact_angle_deg, box_length):
    """Waiting time for the balance pressure.

    Parameters
    ----------
    contact_angle_deg : float
        contact_angle [degrees]
    box_length : float
        box length [nm]
    """
    contact_angle_rad = contact_angle_deg * np.pi / 180
    h = (2 - np.cos(contact_angle_rad)) * np.cos(contact_angle_rad / 2) ** 4

    prefac = (
        kappa_0_2D(contact_angle_deg)
        * 6
        * box_length**2
        * (kappa_0_3D * box_length**3) ** -h
    )

    return prefac ** (1 / (h - 1))


def cavitation_rate_composition(pressure, contact_angle_deg, box_length, A=None):
    """k_tot = k_3D + k_2D"""

    if A is None:
        A = 6 * box_length**2

    V = box_length**3

    k_0_2D = kappa_0_2D(contact_angle_deg) * A
    k_0_3D = kappa_0_spce * V

    contact_angle_rad = contact_angle_deg * np.pi / 180

    current_rate_2D = cavitation_rate(
        pressure=pressure, k_0=k_0_2D, contact_angle=contact_angle_rad
    )

    current_rate_3D = cavitation_rate(pressure=pressure, k_0=k_0_3D, contact_angle=None)

    return current_rate_3D + current_rate_2D


def cavitation_pressure(
    contact_angle_deg,
    waiting_time,
    size_factor,
    surface_tension=55.0,
    temperature=300,
):
    """
    size_factor : float
        prefactor for the cavitation rate. Can either be an area if contact_angle_deg is not None or an
        volume.
    """
    prefac = 16 * np.pi * (surface_tension * scipy.constants.milli) ** 3
    prefac /= 3 * scipy.constants.Boltzmann * temperature

    if contact_angle_deg is not None:
        contact_angle_rad = contact_angle_deg * np.pi / 180
        geometric_factor = (2 - np.cos(contact_angle_rad)) * np.cos(
            contact_angle_rad / 2
        ) ** 4
        log = np.log(size_factor * kappa_0_2D(contact_angle_deg) * waiting_time)
    else:
        geometric_factor = 1
        log = np.log(size_factor * kappa_0_3D * waiting_time)

    return np.sqrt(prefac * geometric_factor / log) / scipy.constants.bar


def cavitation_pressure_composite(
    contact_angle_deg: float,
    waiting_time: float,
    box_length: float,
    A: float = None,
):
    def fun(pressure, contact_angle_deg, box_length, A):
        return (
            cavitation_rate_composition(pressure, contact_angle_deg, box_length, A)
            - 1 / waiting_time
        )

    if A is None:
        A = 6 * box_length**2

    V = box_length**3

    # for inital guesses taken from cavitation pressure of indiviudal models
    x0 = cavitation_pressure(
        contact_angle_deg=contact_angle_deg,
        waiting_time=waiting_time,
        size_factor=A,
    )
    x1 = cavitation_pressure(
        contact_angle_deg=None,
        waiting_time=waiting_time,
        size_factor=V,
    )

    return scipy.optimize.root_scalar(
        f=fun, x0=x0, x1=x1, args=(contact_angle_deg, box_length, A)
    ).root


def critical_box_size(contact_angle_deg, waiting_time):
    """critical_box_size as given waiting time

    Parameters
    ----------
    contact_angle_deg : float
        contact_angle [degrees]
    waiting_time : float
        time [ns]
    """
    contact_angle_rad = contact_angle_deg * np.pi / 180
    h = (2 - np.cos(contact_angle_rad)) * np.cos(contact_angle_rad / 2) ** 4

    prefac = (
        kappa_0_2D(contact_angle_deg)
        * 6
        * waiting_time
        * (kappa_0_3D * waiting_time) ** -h
    )

    return prefac ** (1 / (3 * h - 2))


def critical_box_size_analytical(box_length, waiting_time):
    """critical_box_size as given waiting time

    Parameters
    ----------
    box_length : float
        box length [nm]
    waiting_time : float
        time [ns]

    Returns
    -------
    theta_star : float
        crossover contact angle in [deg]
    """
    chi_rad = chi / 180 * np.pi

    log_3D = np.log(box_length**3 * waiting_time * kappa_0_3D)
    log_2D = np.log(6 * box_length**2 * waiting_time * chi_rad)

    kernel = (log_2D - log_3D) / (theta_0**-4 - 3 / 16 * log_3D)

    return kernel**0.25 * 180 / np.pi


def critical_radius(contact_angle_deg, box_length, waiting_time):
    """critical_box_size as given waiting time

    Parameters
    ----------
    contact_angle_deg : float
        contact_angle [degrees]
    box_length : float
        box length [nm]
    waiting_time : float
        time [ns]
    """
    contact_angle_rad = contact_angle_deg * np.pi / 180
    h = (2 - np.cos(contact_angle_rad)) * np.cos(contact_angle_rad / 2) ** 4

    return np.sqrt(
        (kappa_0_3D * waiting_time * box_length**3) ** h
        / kappa_0_2D(contact_angle_deg)
        / np.pi
        / waiting_time
    )
