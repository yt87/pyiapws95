Overview
========

**pyiapws95** provides two top-level functions: `water_props` and `saturation_props`.
The former returns single-state water or gas properties::

  temperature [K]
  pressure [Pa]
  density [kg m^-3]
  entropy [Joule kg^-1 K^-1] 
  internal_energy [Joule kg^-1]
  enthalpy [Joule kg^-1]
  gibbs_free_energy [Joule kg^-1]
  isochoric_heat_capacity [Joule kg^-1 K^-1]
  isobaric_heat_capacity [Joule kg^-1 K^-1]
  speed_of_sound [m s^-1]

The latter returns a tuple of above properties for liquid and vapour phases.
Input parameters should be in SI units. When `pint` is available, `pint.Quantity`
can be passed as well. Output is always in SI units, either as `pint.Quantity`
values when an input parameter `units` is set (it is set by default), ot plain
floats.

Examples
---------

Compute thermodynamic water properties for pressure 100 000 Pa and temperature
280 deg K.

.. code-block:: python

  >>> import pyiapws95 as wp  # Wagner and Pruss
  >>> wp.water_props(1e5, 280, units=True)

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

Compute saturation properties for temperature 30 deg C:

.. code-block:: python

  >>> import pyiapws95 as wp
  >>> from pint import UnitRegistry
  >>> ur = UnitRegistry()
  >>> t = ur.Quantity(30, ur.celsius)
  >>> wp.saturation_props(t, units=False)
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
  
When speed is an issue, unit processing can be turned off.
Performance gain is substantial:

.. code-block:: python

  >>> import pyiapws95 as wp
  >>> %timeit wp.water_props(1e5, 280, units=True)
  >>> %timeit wp.water_props(1e5, 280, units=False)

  310 µs ± 2.14 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
  28.5 µs ± 160 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)

The times in this example are for environment with numba installed. Without numba,
performance drops significantly:

.. code-block:: python

  >>> import pyiapws95 as wp
  >>> %timeit wp.water_props(1e5, 280, units=False)

  6.37 ms ± 569 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)

  

