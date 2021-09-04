# Water properties from  Wagner und Pruss 2002 paper.
# Inspired by Fluidika.

from typing import Tuple, Union
import numpy as np

from . import njit

# Formulas 6.1 - 6.3
T_c = 647.096  # [K]
p_c = 2.2064e7  # [Pa]
rho_c = 322.0  # [kg m^-3]
R = 4.6151805e2  # [J kg^-1 K^-1]

# Table 6.1
no = np.array(
    [
        0,
        -8.32044648201,
        6.6832105268,
        3.00632,
        0.012436,
        0.97315,
        1.27950,
        0.96956,
        0.24873,
    ]
)
gammao = np.array([1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105])


# Table 6.2
c = np.array(
    [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        3,
        3,
        3,
        3,
        4,
        6,
        6,
        6,
        6,
        0,
        0,
        0,
    ]
)

d = np.array(
    [
        0,
        1,
        1,
        1,
        2,
        2,
        3,
        4,
        1,
        1,
        1,
        2,
        2,
        3,
        4,
        4,
        5,
        7,
        9,
        10,
        11,
        13,
        15,
        1,
        2,
        2,
        2,
        3,
        4,
        4,
        4,
        5,
        6,
        6,
        7,
        9,
        9,
        9,
        9,
        9,
        10,
        10,
        12,
        3,
        4,
        4,
        5,
        14,
        3,
        6,
        6,
        6,
        3,
        3,
        3,
    ]
)

t = np.array(
    [
        0,
        -0.5,
        0.875,
        1,
        0.5,
        0.75,
        0.375,
        1,
        4,
        6,
        12,
        1,
        5,
        4,
        2,
        13,
        9,
        3,
        4,
        11,
        4,
        13,
        1,
        7,
        1,
        9,
        10,
        10,
        3,
        7,
        10,
        10,
        6,
        10,
        10,
        1,
        2,
        3,
        4,
        8,
        6,
        9,
        8,
        16,
        22,
        23,
        23,
        10,
        50,
        44,
        46,
        50,
        0,
        1,
        4,
    ]
)

n = np.array(
    [
        0,
        0.12533547935523e-01,
        0.78957634722828e01,
        -0.87803203303561e01,
        0.31802509345418,
        -0.26145533859358,
        -0.78199751687981e-02,
        0.88089493102134e-02,
        -0.66856572307965,
        0.20433810950965,
        -0.66212605039687e-04,
        -0.19232721156002,
        -0.25709043003438,
        0.16074868486251,
        -0.40092828925807e-01,
        0.39343422603254e-06,
        -0.75941377088144e-05,
        0.56250979351888e-03,
        -0.15608652257135e-04,
        0.11537996422951e-08,
        0.36582165144204e-06,
        -0.13251180074668e-11,
        -0.62639586912454e-09,
        -0.10793600908932,
        0.17611491008752e-01,
        0.22132295167546,
        -0.40247669763528,
        0.58083399985759,
        0.49969146990806e-02,
        -0.31358700712549e-01,
        -0.74315929710341,
        0.47807329915480,
        0.20527940895948e-01,
        -0.13636435110343,
        0.14180634400617e-01,
        0.83326504880713e-02,
        -0.29052336009585e-01,
        0.38615085574206e-01,
        -0.20393486513704e-01,
        -0.16554050063734e-02,
        0.19955571979541e-02,
        0.15870308324157e-03,
        -0.16388568342530e-04,
        0.43613615723811e-01,
        0.34994005463765e-01,
        -0.76788197844621e-01,
        0.22446277332006e-01,
        -0.62689710414685e-04,
        -0.55711118565645e-09,
        -0.19905718354408,
        0.31777497330738,
        -0.11841182425981,
        -0.31306260323435e02,
        0.31546140237781e02,
        -0.25213154341695e04,
        -0.14874640856724,
        0.31806110878444,
    ]
)

alpha = np.array([20, 20, 20])
beta = np.array([150, 150, 250])
gamma = np.array([1.21, 1.21, 1.25])
epsilon = np.array([1, 1, 1])
a = np.array([3.5, 3.5])
b = np.array([0.85, 0.95])
A = np.array([0.32, 0.32])
B = np.array([0.2, 0.2])
C = np.array([28, 32])
F = np.array([700, 800])  # D has been replaced by F to avoid conflicts
E = np.array([0.3, 0.3])

# Helmholtz Free Energy and its derivatives
HFEValues = Tuple[float, float, float, float, float, float]
# Returned water properties
_Props = Tuple[float, float, float, float, float, float, float]
Props = Tuple[float, float, float, float, float, float, float, float, float, float]


# Table 6.4
@njit
def fun_phio(delta: float, tau: float) -> HFEValues:
    phio = np.log(delta) + no[1] + no[2] * tau + no[3] * np.log(tau)
    phio_d = 1.0 / delta
    phio_t = no[2] + no[3] / tau
    phio_dd = -1.0 / (delta * delta)
    phio_dt = 0.0
    phio_tt = -no[3] / (tau * tau)

    for i in range(4, 9):
        j = i - 4
        ee = np.exp(gammao[j] * tau)
        phio += no[i] * np.log(1 - 1.0 / ee)
        phio_t += no[i] * gammao[j] / (ee - 1)
        phio_tt -= no[i] * ee * (gammao[j] / (ee - 1)) ** 2

    phio_dt = 0.0

    return phio, phio_d, phio_t, phio_dd, phio_dt, phio_tt


# Table 6.5
@njit
def fun_phir(delta: float, tau: float) -> HFEValues:
    phir = 0.0
    phir_d = 0.0
    phir_t = 0.0
    phir_dd = 0.0
    phir_dt = 0.0
    phir_tt = 0.0

    for i in range(1, 8):
        _A = n[i] * delta ** d[i] * tau ** t[i]  # underscore to avoid conflict
        A_d = d[i] / delta * _A
        A_t = t[i] / tau * _A
        A_dd = (d[i] - 1) / delta * A_d
        A_dt = t[i] * d[i] / (tau * delta) * _A
        A_tt = (t[i] - 1) / tau * A_t

        phir += _A
        phir_d += A_d
        phir_t += A_t
        phir_dd += A_dd
        phir_dt += A_dt
        phir_tt += A_tt

    for i in range(8, 52):
        dci = delta ** c[i]
        _B = n[i] * delta ** d[i] * tau ** t[i] * np.exp(-dci)
        B_d = (d[i] - c[i] * dci) / delta * _B
        B_t = t[i] / tau * _B
        B_dd = (d[i] - c[i] * dci - 1) / delta * B_d - dci * (c[i] / delta) ** 2 * _B
        B_dt = t[i] / tau * B_d
        B_tt = (t[i] - 1) / tau * B_t

        phir += _B
        phir_d += B_d
        phir_t += B_t
        phir_dd += B_dd
        phir_dt += B_dt
        phir_tt += B_tt

    for i in range(52, 55):
        j = i - 52

        aux1d = d[i] / delta - 2 * alpha[j] * (delta - epsilon[j])
        aux1t = t[i] / tau - 2 * beta[j] * (tau - gamma[j])
        aux2d = d[i] / (delta * delta) + 2 * alpha[j]
        aux2t = t[i] / (tau * tau) + 2 * beta[j]

        _C = (
            n[i]
            * delta ** d[i]
            * tau ** t[i]
            * np.exp(-alpha[j] * (delta - epsilon[j]) ** 2 - beta[j] * (tau - gamma[j]) ** 2)
        )
        C_d = aux1d * _C
        C_t = aux1t * _C
        C_dd = aux1d * C_d - aux2d * _C
        C_dt = aux1d * aux1t * _C
        C_tt = aux1t * C_t - aux2t * _C

        phir += _C
        phir_d += C_d
        phir_t += C_t
        phir_dd += C_dd
        phir_dt += C_dt
        phir_tt += C_tt

    for i in range(55, 57):
        j = i - 55
        dd = (delta - 1) ** 2
        tt = (tau - 1) ** 2

        theta = (1 - tau) + A[j] * dd ** (0.5 / E[j])
        theta_d = (theta + tau - 1) / (delta - 1) / E[j]
        theta_dd = (1.0 / E[j] - 1) * theta_d / (delta - 1)

        psi = np.exp(-C[j] * dd - F[j] * tt)
        psi_d = -2 * C[j] * (delta - 1) * psi
        psi_t = -2 * F[j] * (tau - 1) * psi
        psi_dd = -2 * C[j] * (psi + (delta - 1) * psi_d)
        psi_dt = 4 * C[j] * F[j] * (delta - 1) * (tau - 1) * psi
        psi_tt = -2 * F[j] * (psi + (tau - 1) * psi_t)

        Delta = theta * theta + B[j] * dd ** a[j]
        Delta_d = 2 * (theta * theta_d + a[j] * (Delta - theta * theta) / (delta - 1))
        Delta_t = -2 * theta
        Delta_dd = 2 * (
            theta_d * theta_d
            + theta * theta_dd
            + a[j]
            * (
                (Delta_d - 2 * theta * theta_d) / (delta - 1)
                - (Delta - theta * theta) / (delta - 1) ** 2
            )
        )
        Delta_dt = -2 * theta_d
        Delta_tt = 2

        DeltaPow = Delta ** b[j]
        DeltaPow_d = b[j] * Delta_d / Delta * DeltaPow
        DeltaPow_t = b[j] * Delta_t / Delta * DeltaPow
        DeltaPow_dd = (
            b[j] * Delta_dd / Delta + b[j] * (b[j] - 1) * (Delta_d / Delta) ** 2
        ) * DeltaPow
        DeltaPow_dt = (
            b[j] * Delta_dt / Delta + b[j] * (b[j] - 1) * (Delta_d / Delta) * (Delta_t / Delta)
        ) * DeltaPow
        DeltaPow_tt = (
            b[j] * Delta_tt / Delta + b[j] * (b[j] - 1) * (Delta_t / Delta) ** 2
        ) * DeltaPow

        D = n[i] * DeltaPow * delta * psi
        D_d = n[i] * (DeltaPow * (psi + delta * psi_d) + DeltaPow_d * delta * psi)
        D_t = n[i] * delta * (DeltaPow_t * psi + DeltaPow * psi_t)
        D_dd = n[i] * (
            DeltaPow * (2 * psi_d + delta * psi_dd)
            + 2 * DeltaPow_d * (psi + delta * psi_d)
            + DeltaPow_dd * delta * psi
        )
        D_dt = n[i] * (
            DeltaPow * (psi_t + delta * psi_dt)
            + delta * DeltaPow_d * psi_t
            + DeltaPow_t * (psi + delta * psi_d)
            + DeltaPow_dt * delta * psi
        )
        D_tt = n[i] * delta * (DeltaPow_tt * psi + 2 * DeltaPow_t * psi_t + DeltaPow * psi_tt)

        phir += D
        phir_d += D_d
        phir_t += D_t
        phir_dd += D_dd
        phir_dt += D_dt
        phir_tt += D_tt

    return phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt


# Formula 6.4
@njit
def fun_phi(delta: float, tau: float) -> HFEValues:
    phio, phio_d, phio_t, phio_dd, phio_dt, phio_tt = fun_phio(delta, tau)
    phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt = fun_phir(delta, tau)

    phi = phio + phir
    phi_d = phio_d + phir_d
    phi_t = phio_t + phir_t
    phi_dd = phio_dd + phir_dd
    phi_dt = phio_dt + phir_dt
    phi_tt = phio_tt + phir_tt

    return phi, phi_d, phi_t, phi_dd, phi_dt, phi_tt


# Table 6.3
@njit
def _water_props(delta: float, tau: float) -> _Props:
    T = T_c / tau
    phio, phio_d, phio_t, phio_dd, phio_dt, phio_tt = fun_phio(delta, tau)
    phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt = fun_phir(delta, tau)

    phi = phio + phir
    t_phi_t = (phio_t + phir_t) * tau
    tt_phi_tt = (phio_tt + phir_tt) * tau * tau
    d_phir_d = phir_d * delta
    dd_phir_dd = phir_dd * delta * delta
    dt_phir_dt = phir_dt * delta * tau

    s = R * (t_phi_t - phi)
    u = R * T * t_phi_t
    h = R * T * (1 + t_phi_t + d_phir_d)
    g = h - T * s
    cv = -R * tt_phi_tt

    aux1 = (1 + d_phir_d - dt_phir_dt) ** 2
    aux2 = 1 + 2 * d_phir_d + dd_phir_dd

    cp = cv + R * aux1 / aux2
    w = np.sqrt(R * T * (aux2 - aux1 / tt_phi_tt))

    return s, u, h, g, cv, cp, w


# Find properties given pressure, temperature and an estimate of density
@njit
def water_props(pressure: float, temperature: float, density: float) -> Union[Props, None]:
    """Computes water properties.

    Implements formulas from table 6.3 in
    `Wagner and Pruss <https://aip.scitation.org/doi/10.1063/1.1461829>`_

    Parameters
    ----------
    pressure : float
        Pressure [Pa]
    temperature : float
        Temperature [K]
    density : float
        Initial density [kg/m^3]

    Returns
    -------
    tuple or None
        Calculated values: 'temperature', 'pressure', 'density', 'entropy',
        'internal_energy', 'enthalpy', 'gibbs_free_energy', 'isochoric_heat_capacity',
        'isobaric_heat_capacity', 'speed_of_sound'. When iteration fails to converge
        `None` is returned.
    """
    max_iters = 100
    tolerance = 1e-7 * rho_c

    delta = density / rho_c
    tau = T_c / temperature
    c = rho_c * R * temperature
    # Calculate density from the first equation in Table 6.3.
    for k in range(max_iters):
        phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt = fun_phir(delta, tau)
        d_phir_d = phir_d * delta
        dd_phir_dd = phir_dd * delta * delta
        f = c * delta * (1 + d_phir_d) - pressure
        if abs(f) < tolerance:
            # Found density, calculate remaining properties
            props = (temperature, pressure, delta * rho_c)
            return props + _water_props(delta, tau)
        df = c * (1 + 2 * d_phir_d + dd_phir_dd)
        delta -= f / df

    # Iteration failed
    return None


# Find properties at saturation given temperature
# Formulas 6.9a - 6.9c
@njit
def saturation_props(temperature: float) -> Union[Tuple[Props, Props], None]:
    """Computes water properties at saturation pressure.

    Implements formulas 6.9a - 6.9c from
    `Wagner and Pruss <https://aip.scitation.org/doi/10.1063/1.1461829>`_

    Parameters
    ----------
    temperature : float
        Temperature [K]

    Returns
    -------
    tuple or None
        Tuple with fields: 'liquid' and 'vapour'. Each field is a tuple with water
        properties for the relevant phase. When iteration fails to converge
        `None` is returned.

    See Also
    --------
    water_props : Return water properties for single phase.
    """
    max_iters = 100
    tolerance = 1e-6 * rho_c

    tau = T_c / temperature
    c = rho_c * R * temperature
    delta1 = saturated_liquid_density(temperature) / rho_c
    delta2 = saturated_vapour_density(temperature) / rho_c
    ps = saturated_vapour_pressure(temperature)
    X = np.array([delta1, delta2, ps])
    for k in range(max_iters):
        delta1, delta2, ps = X
        phir1, phir_d1, phir_t1, phir_dd1, phir_dt1, phir_tt1 = fun_phir(delta1, tau)
        phir2, phir_d2, phir_t2, phir_dd2, phir_dt2, phir_tt2 = fun_phir(delta2, tau)

        d_phir_d1 = delta1 * phir_d1
        dd_phir_dd1 = delta1 * delta1 * phir_dd1
        d_phir_d2 = delta2 * phir_d2
        dd_phir_dd2 = delta2 * delta2 * phir_dd2

        f1 = c * delta1 * (1 + d_phir_d1) - ps
        f2 = c * delta2 * (1 + d_phir_d2) - ps
        # By subtraction of 6.9a and 6.9.b from 6.9c
        f3 = (phir1 - phir2 + np.log(delta1 / delta2)) + d_phir_d1 - d_phir_d2

        if max(abs(f1), abs(f2)) < tolerance:
            props = (temperature, ps)
            props_liquid = props + (delta1 * rho_c,) + _water_props(delta1, tau)
            props_vapour = props + (delta2 * rho_c,) + _water_props(delta2, tau)
            return props_liquid, props_vapour

        F = np.array([f1, f2, f3])
        aux1 = 1 + 2 * d_phir_d1 + dd_phir_dd1
        aux2 = 1 + 2 * d_phir_d2 + dd_phir_dd2
        J = np.array([[c * aux1, 0, -1], [0, c * aux2, -1], [aux1 / delta1, -aux2 / delta2, 0]])
        D = np.linalg.solve(J, F)
        X -= D

    return None


# Formula 2.5
@njit
def saturated_vapour_pressure(
    temperature: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """Computes saturation vapour pressure.

    Implements formula 2.5 from
    `Wagner and Pruss <https://aip.scitation.org/doi/10.1063/1.1461829>`_

    Parameters
    ----------
    temperature : float or numpy.ndarray
        Temperature [K]

    Returns
    -------
    float or numpy.ndarray
        Vapour pressure [Pa]
    """

    a1 = -7.85951783
    a2 = 1.84408259
    a3 = -11.7866497
    a4 = 22.6807411
    a5 = -15.9618719
    a6 = 1.80122502

    t = 1.0 - temperature / T_c
    t15 = t ** 1.5
    t30 = t15 * t15
    t35 = t15 * t * t
    t40 = t30 * t
    t75 = t35 * t40

    e = a1 * t + a2 * t15 + a3 * t30 + a4 * t35 + a5 * t40 + a6 * t75
    return p_c * np.exp(T_c / temperature * e)


# Formula 2.6
@njit
def saturated_liquid_density(temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Computes saturated liquid density.

    Implements formula 2.6 from
    `Wagner and Pruss <https://aip.scitation.org/doi/10.1063/1.1461829>`_

    Parameters
    ----------
    temperature : float or numpy.ndarray
        Temperature [K]

    Returns
    -------
    float or numpy.ndarray
        Density [kg m^-3]
    """
    b1 = 1.99274064
    b2 = 1.09965342
    b3 = -0.510839303
    b4 = -1.75493479
    b5 = -45.5170352
    b6 = -6.74694450e05

    t = 1.0 - temperature / T_c
    t13 = t ** (1.0 / 3.0)
    t23 = t13 * t13
    t53 = t13 * t23 * t23
    t163 = t13 * t53 * t53 * t53
    t433 = t163 * t163 * t53 * t * t
    t1103 = t433 * t433 * t163 * t53 * t

    e = 1.0 + b1 * t13 + b2 * t23 + b3 * t53 + b4 * t163 + b5 * t433 + b6 * t1103
    return rho_c * e


# Formula 2.7
@njit
def saturated_vapour_density(temperature: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """Computes saturated vapour density.

    Implements formula 2.7 from
    `Wagner and Pruss <https://aip.scitation.org/doi/10.1063/1.1461829>`_

    Parameters
    ----------
    temperature : float or numpy.ndarray
        Temperature [K]

    Returns
    -------
    float or numpy.ndarray
        Density [kg m^-3]
    """
    c1 = -2.03150240
    c2 = -2.68302940
    c3 = -5.38626492
    c4 = -17.2991605
    c5 = -44.7586581
    c6 = -63.9201063

    t = 1.0 - temperature / T_c
    t16 = t ** (1.0 / 6.0)
    t26 = t16 * t16
    t46 = t26 * t26
    t86 = t46 * t46
    t186 = t86 * t86 * t26
    t376 = t186 * t186 * t16
    t716 = t376 * t186 * t86 * t86

    e = c1 * t26 + c2 * t46 + c3 * t86 + c4 * t186 + c5 * t376 + c6 * t716
    return rho_c * np.exp(e)
