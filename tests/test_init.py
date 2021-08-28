import pint
import pytest

from pyiapws95 import wagner, water_props, saturation_props


def test_saturation_props_liquid():
    tK = 300
    ur = pint.UnitRegistry()
    tC_units = (tK * ur.kelvin).to(ur.celsius)

    propsK = saturation_props(tK, units=False).liquid
    propsC = saturation_props(tC_units, units=True).liquid

    assert propsK.temperature == pytest.approx(propsC.temperature.to("kelvin").m)
    assert propsK.pressure == pytest.approx(propsC.pressure.m)
    assert propsK.density == pytest.approx(propsC.density.m)
    assert propsK.entropy == pytest.approx(propsC.entropy.m)
    assert propsK.enthalpy == pytest.approx(propsC.enthalpy.m)
    assert propsK.isobaric_heat_capacity == pytest.approx(propsC.isobaric_heat_capacity.m)
    assert propsK.isochoric_heat_capacity == pytest.approx(propsC.isochoric_heat_capacity.m)
    assert propsK.speed_of_sound == pytest.approx(propsC.speed_of_sound.m)


def test_saturation_props_vapour():
    tK = 500
    ur = pint.UnitRegistry()
    tC_units = (tK * ur.kelvin).to(ur.celsius)

    propsK = saturation_props(tK, units=False).vapour
    propsC = saturation_props(tC_units, units=True).vapour

    assert propsK.temperature == pytest.approx(propsC.temperature.to("kelvin").m)
    assert propsK.pressure == pytest.approx(propsC.pressure.m)
    assert propsK.density == pytest.approx(propsC.density.m)
    assert propsK.entropy == pytest.approx(propsC.entropy.m)
    assert propsK.enthalpy == pytest.approx(propsC.enthalpy.m)
    assert propsK.isobaric_heat_capacity == pytest.approx(propsC.isobaric_heat_capacity.m)
    assert propsK.isochoric_heat_capacity == pytest.approx(propsC.isochoric_heat_capacity.m)
    assert propsK.speed_of_sound == pytest.approx(propsC.speed_of_sound.m)


def test_water_props():
    pPa = 1e06
    tK = 500
    ur = pint.UnitRegistry()
    pMPa_units = (pPa * ur.pascal).to(ur.megapascal)
    tC_units = (tK * ur.kelvin).to(ur.celsius)

    propsK = water_props(pPa, tK, units=False)
    propsC = water_props(pMPa_units, tC_units, units=True)

    assert propsK.temperature == pytest.approx(propsC.temperature.to("kelvin").m)
    assert propsK.pressure == pytest.approx(propsC.pressure.to("pascal").m)
    assert propsK.density == pytest.approx(propsC.density.m)
    assert propsK.entropy == pytest.approx(propsC.entropy.m)
    assert propsK.enthalpy == pytest.approx(propsC.enthalpy.m)
    assert propsK.isobaric_heat_capacity == pytest.approx(propsC.isobaric_heat_capacity.m)
