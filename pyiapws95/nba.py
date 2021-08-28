# Numba implementation

import numpy as np
from numba import float32, float64, njit, vectorize

R_star = 8.31432e03
R = 4.6151805e02
M = 2.89655e01
EPS = R_star / M / R


@njit([float32(float32), float64(float64)])
def enthalpy_dry_air(t):
    M = 2.89655e01
    A = 2.93959e01
    B = -2.76362e00
    C = 4.57385e00
    D = 5.34084e00
    E = -1.45350e-03
    F = -7.97022e00  # adjusted to h(273.15) = 0

    x = t / 1000.0
    h = 1e06 * ((((D / 4 * x + C / 3) * x + B / 2) * x + A) * x - E / x + F) / M
    dh_dt = 1e03 * (((D * x + C) * x + B) * x + A + E / (x * x)) / M

    return h, dh_dt


@njit([float32(float32), float64(float64)])
def enthalpy_saturated_vapour(t):
    a = -2.6630873e-05
    b = -2.9011536e-03
    c = -3.2596599e-01
    d = 1.8347209e03
    e = 2.50089686e06

    x = t - 273.15
    h = (((a * x + b) * x + c) * x + d) * x + e
    dh_dt = ((4 * a * x + 3 * b) * x + 2 * c) * x + d

    return h, dh_dt


@njit([float32(float32, float32), float64(float64, float64)])
def enthalpy_water(t, p):
    A = 2.62e02
    B = 4.5e-01
    a = 1.15868908e-03
    b = -4.64884789e-05
    c = 9.63993632e-02
    d = -3.36641920e00
    e = 1.83541646e02
    f = 5.69684155e03
    g = -2.90951343e03
    h = -4.58629080e04

    x = (t - A) ** B
    y = p
    h = a * y + b * x * y + ((((c * x + d) * x + e) * x + f) * x + g) * x + h
    dh_dt = (b * y + (((5 * c * x + 4 * d) * x + 3 * e) * x + 2 * f) * x + g) * B * x / (t - A)

    return h, dh_dt


@njit([float32(float32), float64(float64)])
def saturation_vapour_pressure(t):
    a = 3.4494e01
    b = 4.92499e03
    c = 1.57
    d1 = 3.605e01
    d2 = 1.6815e02
    es = np.exp(a - b / (t - d1)) / (t - d2) ** c
    des_dt = es * (b / (t - d1) ** 2 - c / (t - d2))
    return es, des_dt


@njit([float32(float32, float32), float64(float64, float64)])
def saturation_mixing_ratio(t, p):
    es, des = saturation_vapour_pressure(t)
    ws = EPS * es / (p - es)
    dws = EPS * des * p / (p - es) ** 2
    return ws, dws


@njit(
    [
        float32(float32, float32, float32, float32),
        float64(float64, float64, float64, float64),
    ]
)
def newton(t, p, tdb, w):
    h, dh = enthalpy_dry_air(tdb)
    hg, dhg = enthalpy_saturated_vapour(tdb)
    hw, dhw = enthalpy_water(t, p)
    ws, dws = saturation_mixing_ratio(t, p)
    h1, dh1 = enthalpy_dry_air(t)
    hg1, dhg1 = enthalpy_saturated_vapour(t)
    f = (h1 + ws * hg1 - (h + w * hg + (ws - w) * hw),)
    df_dt = dh1 + dws * hg1 + ws * dhg1 - (dws * hw + (ws - w) * dhw)
    return f / df_dt


@vectorize(
    [
        float32(float32, float32, float32, float32),
        float64(float64, float64, float64, float64),
    ],
    target="parallel",
)
def wet_bulb_temperature(temperature, relative_humidity, pressure, precision=1e-2):
    """
    Computes isobaric wet bulb temperature.

    Parameters
    ----------
    temperature : array_like
        Dry bulb temperature [deg K]
    relative_humidity : array_like
        Relative humidity [fraction]
    pressure : array_like
        Air pressure [Pa]
    precision : float
        Maximum absolute error [deg K]

    Returns
    -------
    float or ndarray
        Wet bulb temperature [deg K]
    """
    niter = 10
    eps = 4 * np.sqrt(precision)
    es, des_dt = saturation_vapour_pressure(temperature)
    e = es * relative_humidity
    w = EPS * e / (pressure - e)
    wet_bulb_temperature = temperature
    for i in range(niter):
        delta = newton(wet_bulb_temperature, pressure, temperature, w)
        wet_bulb_temperature -= delta
        if abs(delta) < eps:
            break
    return wet_bulb_temperature
