Overview
========

**pyiapws95** provides two top-level functions: :py:func:`~pyiapws95.water_props`
and :py:func:`~pyiapws95.saturation_props`.
The former returns single-state water or gas properties::

.. parsed-literal::

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
Input parameters should be in SI units. When pint is available, pint.Quantity
can be passed as well. Output is always in SI units, either as pint.Quantity
values when an input parameter units is set (it is set by default), ot plain
floats.

Examples
---------

Compute thermodynamic water properties for pressure 100 000 Pa and temperature
280 deg K.

.. code-block:: python

   >>> import pyiapws95 as wp  # Wagner and Pruss or Water Properties
   >>> wp.water_props(1e5, 280, units=True)

.. parsed-literal::

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
   >>> wp.saturation_props(t, units=True)

.. parsed-literal::

   Liquid WaterProps(
     temperature = 303.15 kelvin
     pressure = 4246.970836833437 pascal
     density = 995.6061774738458 kilogram / meter ** 3
     entropy = 436.7548587214662 joule / kelvin / kilogram
     internal_energy = 125729.70746867526 joule / kilogram
     enthalpy = 125733.97318230748 joule / kilogram
     gibbs_free_energy = -6668.262239104966 joule / kilogram
     isochoric_heat_capacity = 4117.534162924182 joule / kelvin / kilogram
     isobaric_heat_capacity = 4180.08366959557 joule / kelvin / kilogram
     speed_of_sound = 1508.9888796749242 meter / second
   )
   Vapour WaterProps(
     temperature = 303.15 kelvin
     pressure = 4246.970836833437 pascal
     density = 0.030415211808981893 kilogram / meter ** 3
     entropy = 8451.965908946451 joule / kelvin / kilogram
     internal_energy = 2415912.08580385 joule / kilogram
     enthalpy = 2555545.203058002 joule / kilogram
     gibbs_free_energy = -6668.262239113916 joule / kilogram
     isochoric_heat_capacity = 1445.240733936132 joule / kelvin / kilogram
     isobaric_heat_capacity = 1918.0102698516084 joule / kelvin / kilogram
     speed_of_sound = 430.0307654967225 meter / second
   )
  
When speed is an issue, unit processing can be turned off.
Performance gain is substantial:

.. code-block:: python

   >>> import pyiapws95 as wp
   >>> %timeit wp.water_props(1e5, 280, units=True)
   >>> %timeit wp.water_props(1e5, 280, units=False)

.. parsed-literal::

   310 µs ± 2.14 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)
   28.5 µs ± 160 ns per loop (mean ± std. dev. of 7 runs, 10000 loops each)

The times in this example are for an environment with numba installed. Without numba,
performance drops significantly:

.. code-block:: python

   >>> import pyiapws95 as wp
   >>> %timeit wp.water_props(1e5, 280, units=False)

.. parsed-literal::

   6.37 ms ± 569 µs per loop (mean ± std. dev. of 7 runs, 1 loop each)
