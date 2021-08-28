from typing import Iterable, List, NamedTuple, Union

try:
    from pint import UnitRegistry, Quantity

    _r = UnitRegistry()
    _units: Union[None, List["Quantity"]] = [
        _r.kelvin,
        _r.pascal,
        _r.kg / _r.m ** 3,
        _r.joule / _r.kg / _r.kelvin,
        _r.joule / _r.kg,
        _r.joule / _r.kg,
        _r.joule / _r.kg,
        _r.joule / _r.kg / _r.kelvin,
        _r.joule / _r.kg / _r.kelvin,
        _r.m / _r.s,
    ]
except ImportError:
    _units = None

from . import _version
__version__ = _version.get_versions()['version']

from . import wagner, tables

_value_type = Union["Quantity", float]


class WaterProps(NamedTuple):
    temperature: _value_type
    pressure: _value_type
    density: _value_type
    entropy: _value_type
    internal_energy: _value_type
    enthalpy: _value_type
    gibbs_free_energy: _value_type
    isochoric_heat_capacity: _value_type
    isobaric_heat_capacity: _value_type
    speed_of_sound: _value_type

    def __repr__(self):
        s = "\n".join(["  {:s} = {}".format(n, getattr(self, n)) for n in self._fields])
        return "WaterProps(\n{}\n)".format(s)


class SaturatedProps(NamedTuple):
    liquid: WaterProps
    vapour: WaterProps

    def __repr__(self):
        return "\n".join(["Liquid " + repr(self.liquid), "Vapour " + repr(self.vapour)])


def with_units(props: Iterable[float]) -> List["Quantity"]:
    if not _units:
        raise ValueError("Units not defined")
    return [v * u for (v, u) in zip(props, _units)]


def water_props(pressure, temperature, density=None, units=True):
    """Computes water properties.

    Parameters
    ----------
    pressure : float or pint.Quantity
        Pressure, either in pascals or as a pint Quantity.
    temperature : float or pint.Quantity
        Temperature, in degrees Kelvin or as a pint Quantity.
    density : float, pint.Quantity or None
        Initial density. If None, the density will be determined by lookup
        of Table 13.2. If numeric, must be kg/m^3.
    units : bool
        If pint is avaliable and `units` is True, return pint Quantites.
        Default is True.

    Returns
    -------
    typing.NamedTuple
        Tuple with fields: 'temperature', 'pressure', 'density', 'entropy',
        'internal_energy', 'enthalpy', 'gibbs_free_energy', 'isochoric_heat_capacity',
        'isobaric_heat_capacity', 'speed_of_sound'

    Raises
    ------
    ValueError
        When the algorithm fails to converge.

    See Also
    --------
    saturation_props : Return water properties at saturation.
        
    Examples
    --------
    >>> water_props(1e5, 280)
    WaterProps(
      temperature = 280.0 kelvin
      pressure = 100000.0 pascal
      density = 999.9103569375845 kilogram / meter ** 3
      entropy = 104.11392202805507 joule / kelvin / kilogram
      internal_energy = 28794.08780871546 joule / kilogram
      enthalpy = 28894.09677381184 joule / kilogram
      gibbs_free_energy = -257.8013940435776 joule / kilogram
      isochoric_heat_capacity = 4199.837682609676 joule / kelvin / kilogram
      isobaric_heat_capacity = 4200.944742530026 joule / kelvin / kilogram
      speed_of_sound = 1434.2746295730974 meter / second
    )
    """
    if isinstance(pressure, Quantity):
        pressure = pressure.to("pascal").m
    pressure = float(pressure)
    if isinstance(temperature, Quantity):
        temperature = temperature.to("kelvin").m
    temperature = float(temperature)
    if density is None:
        density = tables.initial_rho(pressure, temperature)
    else:
        if isinstance(density, Quantity):
            density = density.to("kg/m^3").m
        density = float(density)

    props = wagner.water_props(pressure, temperature, density)
    if props is None:
        raise ValueError(
            "Computation failed for p={:e}, t={:e}, rho={:e}".format(
                pressure, temperature, density
            )
        )
    if not units or _units is None:
        return WaterProps(*props)
    else:
        return WaterProps(*with_units(props))


def saturation_props(temperature, units=True):
    """Computes water properties at saturation pressure.

    Parameters
    ----------
    temperature : float or pint.Quantity
        Temperature, in degrees Kelvin or as a pint Quantity.
    units : bool
        If pint is avaliable and `units` is True, return pint Quantites.
        Default is True.

    Returns
    -------
    typing.NamedTuple
        Tuple with fields: 'liquid' and 'vapour'. Each field is a tuple with water
        properties for the relevant phase.

    Raises
    ------
    ValueError
        When the algorithm fails to converge.

    See Also
    --------
    water_props : Return water properties for single phase.
        
    Examples
    --------
    >>> saturation_props(280, units=False)
    Liquid WaterProps(
      temperature = 280.0
      pressure = 991.8203199195548
      density = 999.862210022028
      entropy = 104.11824993199234
      internal_energy = 28795.297190165784
      enthalpy = 28796.28914716716
      gibbs_free_energy = -356.82083379069445
      isochoric_heat_capacity = 4200.257591191563
      isobaric_heat_capacity = 4201.350463857432
      speed_of_sound = 1434.1151160257625
    )
    Vapour WaterProps(
      temperature = 280.0
      pressure = 991.8203199195548
      density = 0.007681162498197408
      entropy = 8977.875870567612
      internal_energy = 2384324.697167183
      enthalpy = 2513448.422925141
      gibbs_free_energy = -356.820833790116
      isochoric_heat_capacity = 1424.168493620011
      isobaric_heat_capacity = 1891.341041100502
      speed_of_sound = 413.9245165642786
    )
    """
    if isinstance(temperature, Quantity):
        temperature = temperature.to("kelvin").m
    temperature = float(temperature)
    props = wagner.saturation_props(temperature)
    if props is None:
        raise ValueError("Computation failed for t={:e} K".format(temperature))
    if not units or _units is None:
        return SaturatedProps(WaterProps(*props[0]), WaterProps(*props[1]))
    else:
        return SaturatedProps(
            WaterProps(*with_units(props[0])), WaterProps(*with_units(props[1]))
        )
