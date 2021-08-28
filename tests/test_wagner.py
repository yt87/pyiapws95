import pytest

from pyiapws95 import wagner, tables


def test_phio500():
    # Table 6.6
    expected_phio = 0.204797743e1
    expected_phio_d = 0.384236747
    expected_phio_t = 0.904611106e1
    expected_phio_dd = -0.147637878
    expected_phio_dt = 0.0
    expected_phio_tt = -0.193249185e1

    T = 500.0
    rho = 838.025

    tau = wagner.T_c / T
    delta = rho / wagner.rho_c

    phio, phio_d, phio_t, phio_dd, phio_dt, phio_tt = wagner.fun_phio(delta, tau)

    assert phio == pytest.approx(expected_phio)
    assert phio_d == pytest.approx(expected_phio_d)
    assert phio_t == pytest.approx(expected_phio_t)
    assert phio_dd == pytest.approx(expected_phio_dd)
    assert phio_dt == pytest.approx(expected_phio_dt)
    assert phio_tt == pytest.approx(expected_phio_tt)


def test_phio647():
    # Table 6.6
    expected_phio = -0.156319605e1
    expected_phio_d = 0.899441341
    expected_phio_t = 0.980343918e1
    expected_phio_dd = -0.808994726
    expected_phio_dt = 0.0
    expected_phio_tt = -0.343316334e1

    T = 647.0
    rho = 358.0

    tau = wagner.T_c / T
    delta = rho / wagner.rho_c

    phio, phio_d, phio_t, phio_dd, phio_dt, phio_tt = wagner.fun_phio(delta, tau)

    assert phio == pytest.approx(expected_phio)
    assert phio_d == pytest.approx(expected_phio_d)
    assert phio_t == pytest.approx(expected_phio_t)
    assert phio_dd == pytest.approx(expected_phio_dd)
    assert phio_dt == pytest.approx(expected_phio_dt)
    assert phio_tt == pytest.approx(expected_phio_tt)


def test_phir500():
    # Table 6.6
    expected_phir = -0.342693206e1
    expected_phir_d = -0.364366650
    expected_phir_t = -0.581403435e1
    expected_phir_dd = 0.856063701
    expected_phir_dt = -0.112176915e1
    expected_phir_tt = -0.223440737e1

    T = 500.0
    rho = 838.025

    tau = wagner.T_c / T
    delta = rho / wagner.rho_c

    phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt = wagner.fun_phir(delta, tau)

    assert phir == pytest.approx(expected_phir)
    assert phir_d == pytest.approx(expected_phir_d)
    assert phir_t == pytest.approx(expected_phir_t)
    assert phir_dd == pytest.approx(expected_phir_dd)
    assert phir_dt == pytest.approx(expected_phir_dt)
    assert phir_tt == pytest.approx(expected_phir_tt)


def test_phir647():
    # Table 6.6
    expected_phir = -0.121202657e1
    expected_phir_d = -0.714012024
    expected_phir_t = -0.321722501e1
    expected_phir_dd = 0.475730696
    expected_phir_dt = -0.133214720e1
    expected_phir_tt = -0.996029507e1

    T = 647.0
    rho = 358.0

    tau = wagner.T_c / T
    delta = rho / wagner.rho_c

    phir, phir_d, phir_t, phir_dd, phir_dt, phir_tt = wagner.fun_phir(delta, tau)

    assert phir == pytest.approx(expected_phir)
    assert phir_d == pytest.approx(expected_phir_d)
    assert phir_t == pytest.approx(expected_phir_t)
    assert phir_dd == pytest.approx(expected_phir_dd)
    assert phir_dt == pytest.approx(expected_phir_dt)
    assert phir_tt == pytest.approx(expected_phir_tt)


def test_saturated_vapour_pressure():
    expected_p = 611.657

    T = 273.16

    p = wagner.saturated_vapour_pressure(T)

    assert p == pytest.approx(expected_p)


def test_saturated_liquid_density():
    expected_rho = 958.365

    T = 373.1243

    rho = wagner.saturated_liquid_density(T)

    assert rho == pytest.approx(expected_rho)


def test_saturated_vapour_density():
    expected_rho = 0.597586

    T = 373.1243

    rho = wagner.saturated_vapour_density(T)

    assert rho == pytest.approx(expected_rho)


def test_water_props():
    rows = [0, 17, 19, 153]

    for i in rows:
        T, p, rho, u, h, s, cv, cp, w = tables.Table13_2[i, :]
        rho_init = tables.initial_rho(p, T)
        props = wagner.water_props(p, T, rho_init)

        assert props[2] == pytest.approx(rho, abs=1e-3)
        assert props[3] == pytest.approx(s, abs=1e-1)
        # FIXME: check this, near boiling
        assert props[4] == pytest.approx(u, abs=2)
        assert props[5] == pytest.approx(h, abs=2)
        assert props[7] == pytest.approx(cv, abs=1e-1)
        assert props[8] == pytest.approx(cp, abs=1e-1)
        assert props[9] == pytest.approx(w, abs=1e-1)


def test_saturation_liquid_props():
    rows = [0, 80, 185]

    for i in rows:
        T, p, rho, h, s, cv, cp, w = tables.Table13_1_liquid[i, :]
        props, props_v = wagner.saturation_props(T)
        assert props[2] == pytest.approx(rho, abs=1e-2)
        assert props[3] == pytest.approx(s, abs=1e-1)
        assert props[5] == pytest.approx(h, abs=2)
        assert props[7] == pytest.approx(cv, rel=1e-4)
        assert props[8] == pytest.approx(cp, rel=1e-4)
        assert props[9] == pytest.approx(w, rel=1e-4)


def test_saturation_vapour_props():
    rows = [0, 80, 185]

    for i in rows:
        T, p, rho, h, s, cv, cp, w = tables.Table13_1_vapour[i, :]
        props_l, props = wagner.saturation_props(T)
        assert props[2] == pytest.approx(rho, abs=1e-2)
        assert props[3] == pytest.approx(s, abs=1e-1)
        assert props[5] == pytest.approx(h, abs=5)
        assert props[7] == pytest.approx(cv, rel=1e-4)
        assert props[8] == pytest.approx(cp, rel=1e-4)
        assert props[9] == pytest.approx(w, rel=1e-4)
