"""Functions to fit the attempt frequency surface density κ_2D.

Resulting function for κ_2D is given in constants.py as kappa_0_2D
"""

import numpy as np
import scipy
from constants import gamma_spce


def surface_tension_pressure(Px, Py, Pz, Lz, DPx=0, DPy=0, DPz=0):
    """Calculates the surface tension [mN/m] for an interface in the xy-plane.
    See J. Chem. Phys. 126, 154707 (2007) for details.
    - Px, Py, Py in bar
    - Lz in nm"""
    g = Lz / 2 * (Pz - 0.5 * (Px + Py)) * 1e-1

    dgdPx = Lz / 2
    dgdPy = Lz / 4
    dgdPz = Lz / 4

    Dg = np.sqrt((dgdPx * DPx) ** 2 + (dgdPy * DPy) ** 2 + (dgdPz * DPz) ** 2) * 1e-1

    if Dg != 0:
        return g, Dg

    return g


def contact_angle_energy(contact_angle):
    """Energy contrubution based from the contact angle

    Parameters
    ----------
    contact_angle : float
        contact_angle [rad]
    """

    return (2 - np.cos(contact_angle)) * np.cos(contact_angle / 2) ** 4


def Delta_G(pressure, surface_tension, temperature=300, contact_angle=None):
    """Returns the free energy BARRIER Delta_G^* in k_BT of creating a spherical bubble.

    Parameters
    ----------
    pressure : float or array_like
        pressure [bar]
    surface_tension : float
        surface_tension γ [mN/m]
    temperature : float, optional
        temperature T [K]
    contact_angle : float, optional
        contact_angle [rad]
    """
    energy = 16 * np.pi * (surface_tension * scipy.constants.milli) ** 3
    energy /= (
        3
        * scipy.constants.Boltzmann
        * temperature
        * (pressure * scipy.constants.bar) ** 2
    )

    # Correction for heterogeneous cavitation
    if contact_angle is not None:
        energy *= contact_angle_energy(contact_angle)

    return energy


def cavitation_rate(
    pressure, k_0, surface_tension=gamma_spce, temperature=300, contact_angle=None
):
    """
    Full cavitation rate k either for 3D (contact_angle=None) or 2D.

    Parameters
    ----------
    pressure : float
        in bar
    k_0 : float
        kinetic prefactor in units of 1/time
    surface_tension : float
        surface_tension γ [mN/m]
    temperature : float, optional
        temperature T [K]
    contact_angle : float
        in radians, If none cavitation is 3D
    """
    return k_0 * np.exp(
        -Delta_G(
            pressure,
            surface_tension=surface_tension,
            temperature=temperature,
            contact_angle=contact_angle,
        )
    )


def Delta_G_1D_star(
    pressure, base_radius, surface_tension, contact_angle_defect, temperature=300
):
    """Returns the free energy BARRIER of a pinned bubble.

    Parameters
    ----------
    pressure : float or array_like
        pressure [bar]
    base_radius : float or array_like
        base radius of the bubble on the surface [nm]
    surface_tension : float
        surface_tension γ [mN/m]
    contact_angle_defect : float, optional
        contact_angle of the defect [rad]
    temperature : float, optional
        temperature T [K]
    """

    K = np.sqrt(
        1
        - (
            pressure
            * scipy.constants.bar
            * base_radius
            * scipy.constants.nano
            / (2 * surface_tension * scipy.constants.milli)
        )
        ** 2
    )
    energy = (
        8
        * (surface_tension * scipy.constants.milli) ** 3
        / (pressure * scipy.constants.bar) ** 2
        * (1 + K)
    )

    energy -= (
        surface_tension
        * scipy.constants.milli
        * (base_radius * scipy.constants.nano) ** 2
        * (2 * K - 3 * np.cos(contact_angle_defect))
    )

    return np.pi / 3 * energy / (scipy.constants.Boltzmann * temperature)


def Delta_G_2d(
    base_radius,
    pressure,
    contact_angle,
    surface_tension,
    temperature=300,
):
    """Returns the free energy Delta_G in k_BT of creating a spherical bubble on a surface

    Parameters
    ----------
    base_radius : float or array_like
        base radius of the bubble on the surface [nm]
    pressure : float
        pressure [bar]
    contact_angle : float
        contact_angle of the defect [rad]
    surface_tension : float
        vapor-liquid surface_tension γ [mN/m]
    temperature : float, optional
        temperature T [K]
    """

    # Compute trigometric quantities
    base_area = np.pi * (base_radius * scipy.constants.nano) ** 2
    phi = np.pi - contact_angle
    bubble_radius = base_radius * scipy.constants.nano / np.sin(phi)

    cap_area = 2 * np.pi * bubble_radius**2 * (1 - np.cos(phi))
    bubble_volume = (
        (np.pi / 3) * bubble_radius**3 * (2 - 3 * np.cos(phi) + np.cos(phi) ** 3)
    )

    energy = (
        cap_area * surface_tension * scipy.constants.milli
        + base_area * surface_tension * scipy.constants.milli * np.cos(contact_angle)
        + pressure * scipy.constants.bar * bubble_volume
    )

    return energy / (scipy.constants.Boltzmann * temperature)


def Delta_G_1d(
    bubble_height,
    base_radius,
    pressure,
    contact_angle_defect,
    surface_tension,
    temperature=300,
):
    """Returns the free energy Delta_G in k_BT of creating a spherical bubble on a surface with a defect

    Parameters
    ----------
    bubble_height : float or array_like
        height of the bubble on the defect [nm]
    base_radius : float or array_like
        base radius of the bubble on the surface [nm]
    pressure : float
        pressure [bar]
    contact_angle : float
        contact_angle [rad]
    surface_tension : float
        vapor-liquid surface_tension γ [mN/m]
    temperature : float, optional
        temperature T [K]
    """

    # Compute trigometric quantities
    cos_contact_angle = np.cos(contact_angle_defect)

    base_area = np.pi * (base_radius * scipy.constants.nano) ** 2
    cap_area = np.pi * (
        (base_radius * scipy.constants.nano) ** 2
        + (bubble_height * scipy.constants.nano) ** 2
    )
    bubble_volume = (
        np.pi
        / 6
        * bubble_height
        * scipy.constants.nano
        * (
            3 * (base_radius * scipy.constants.nano) ** 2
            + (bubble_height * scipy.constants.nano) ** 2
        )
    )

    energy = (
        cap_area * surface_tension * scipy.constants.milli
        + base_area * surface_tension * scipy.constants.milli * cos_contact_angle
        + pressure * scipy.constants.bar * bubble_volume
    )

    return energy / (scipy.constants.Boltzmann * temperature)


def Delta_G_lip(pressure, prop_const, adhesion_energy_density, temperature=300):
    """Returns the free energy BARRIER Delta_G in k_BT of creating a a cavity inside a lipid.

    Parameters
    ----------
    pressure : float or array_like
        pressure [bar]
    prop_const : float
        proportionality constant between cavity volume and cavity crosssection α_lip
    temperature : float, optional
        temperature T [K]
    adhesion_energy_density : float
        adhesion energy density ω_lip [kJ/mol/nm^2]
    """
    energy = (
        4
        * (
            adhesion_energy_density
            * scipy.constants.kilo
            / scipy.constants.Avogadro
            / scipy.constants.nano**2
        )
        ** 3
    )
    energy /= (
        27
        * prop_const**2
        * temperature
        * scipy.constants.Boltzmann
        * (pressure * scipy.constants.bar) ** 2
    )

    return energy


def constant_rate_kernel(
    t, rate, k_0, surface_tension, temperature=300, contact_angle=None
):
    """Integral kernel of the constant-rate fit."""

    # use rate instead of pressure here since equations are the same!
    tau_sq = Delta_G(rate, surface_tension, temperature, contact_angle)

    I = t * np.exp(-(tau_sq / t**2))
    I -= np.sqrt(math.pi * tau_sq) * scipy.special.erfc(np.sqrt(tau_sq) / t)
    return np.exp(-k_0 * I)


def p_cav_star(rate, k_0, surface_tension, temperature=300, contact_angle=None, b=1e15):
    """
    Returns $p_{cav}^*(\dot{p})$ as function of the rate by solving the integral kernel.
    Intgration limits are obatained by a bisection algorithm where `b` is right boundary
    for starting the bisection.
    """

    if type(rate) in [list, tuple, np.ndarray]:
        res = np.zeros(len(rate))
        rates = rate
    else:
        res = np.zeros(1)
        rates = np.array([rate])

    res2 = np.zeros(res.shape)
    for i, r in enumerate(rates):
        breakpoint = scipy.optimize.bisect(
            f=lambda t: constant_rate_kernel(
                t, r, k_0, surface_tension, temperature, contact_angle
            )
            - 1e-10,
            a=1e-10,
            b=b,
        )
        t = np.linspace(1e-10, 2 * breakpoint, int(1e4))
        kernel = constant_rate_kernel(
            t, r, k_0, surface_tension, temperature, contact_angle
        )

        res[i] = r * scipy.integrate.simps(kernel, t)

    #     assert_allclose(res, res2)

    if type(rate) in [list, tuple, np.ndarray]:
        return res
    else:
        return res[0]


def fit_kinetiv_prefactor(
    rates,
    p_cav,
    p_cav_err=None,
    surface_tension=None,
    temperature=300,
    contact_angle=None,
    p0=None,
    bounds=None,
):
    """
    Calculates the kinetic prefactor $k_0$ by fitting the mean cavitation pressure $p_{cav}^*$ to
    $$
        p_{cav}^*(\dot{p}) = \int_0^\infty e^{-k_0 I(t)}
    $$
    where $\dot{p}$ is the pressure rate and $I(t)$ the integral kernel.
    If the `surface_tension` $γ$ is not given it is fitted as well.
    For heterogeneous cavitation also the contact angle can be given.

    Parameters
    ----------
    rates : 1-d array
        pressure rates [bar/ns]
    p_cav : 1-d array
        mean cavitation pressure [bar]
    p_cav_err : 1-d array, optional
        error of the mean cavitation pressure [bar]
    surface_tension : float, optional
        surface_tension γ [mN/m]
    temperature : float, optional
        temperature T [K]
    contact_angle : float, optional
        contact_angle [rad]
    p0 : 1-d array , optional
        Initial guess for the parameters.
        If None, then the initial values will all be 1
    bounds : bounds2-tuple of array_like, optional
        Lower and upper bounds on parameters.
        Default for 2 parameterfit [0, 1e5] and for 1 parameter [0, 1e40]

    Returns
    -------
    popt : float
        kinetic prefactor [ns^-1] and obtained surface tension. Present only if `surface_tension = None`.
    p_err : float
        Error of popt [ns^-1]
    """

    if rates.shape[0] != p_cav.shape[0]:
        raise TypeError("expected `rates` and `p_cav` to have same length")

    if p_cav_err is not None:
        if p_cav_err.shape[0] != p_cav.shape[0]:
            raise TypeError("expected `p_cav_err` and `p_cav` to have same length")
        err = 2 / p_cav_err**3
    else:
        err = np.ones(p_cav.shape)

    if surface_tension is None:  # Fit k_0 and gamma
        if bounds is None:
            bounds = (0, 1e5)
        popt, pcov = scipy.optimize.curve_fit(
            f=lambda rate, k_0, surface_tension: 1
            / p_cav_star(np.exp(rate), k_0, surface_tension, temperature, contact_angle)
            ** 2,
            xdata=np.log(rates),
            ydata=1 / p_cav**2,
            sigma=err,
            p0=p0,
            bounds=bounds,
        )
    else:  # One paramater fit
        if bounds is None:
            bounds = (0, 1e40)
        popt, pcov = scipy.optimize.curve_fit(
            f=lambda rate, k_0: 1
            / p_cav_star(np.exp(rate), k_0, surface_tension, temperature, contact_angle)
            ** 2,
            xdata=np.log(rates),
            ydata=1 / p_cav**2,
            sigma=err,
            p0=p0,
            bounds=bounds,
        )

    perr = np.sqrt(np.diag(pcov))

    return popt, perr
