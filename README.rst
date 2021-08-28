pyiapws95
=========

pyiapws95 is a Python implementation of the
`IAPWS formulation 1995  <https://aip.scitation.org/doi/10.1063/1.1461829>`_
for the thermodynamic properties of water. The code is based on a C++ library
`Fluidika <https://phoenix.yizimg.com/reaktoro/fluidika>`_.

The package provides two top-level functions: *water_props* and *saturation_props*.
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

Quick Start
-----------

The only required dependency is numpy. Substantially better performance one can
obtain by installing numba. Pint simplifies unit handling.

The easiest way to try pyiapws95 is to install it with conda in a separate
environment::

  $ conda create -n pyiapws95 -c yt87 python=3.9 numpy numba pint pyiapws95

Then one can obtain thermodynamic properties of water at the sea level and room
temperature as follows:

.. code-block:: python

  >>> from pint import UnitRegistry
  >>> import pyiapws95 as wp
  >>> ur = UnitRegistry()
  >>> pressure = 1013.25 * ur.kilopascals
  >>> temperature = ur.Quantity(18, ur.celsius)
  >>> wp.water_props(pressure, temperature)
  WaterProps(
    temperature = 291.15 kelvin
    pressure = 1013250.0 pascal
    density = 999.0189673750012 kilogram / meter ** 3
    entropy = 267.64473499097977 joule / kelvin / kilogram
    internal_energy = 75487.06425878826 joule / kilogram
    enthalpy = 76501.3092662225 joule / kilogram
    gibbs_free_energy = -1423.4553264012648 joule / kilogram
    isochoric_heat_capacity = 4160.620288359394 joule / kelvin / kilogram
    isobaric_heat_capacity = 4182.645390095095 joule / kelvin / kilogram
    speed_of_sound = 1477.5390482281166 meter / second
  )

Documentation
-------------
HTML documentation is available at https://yt87.github.io/pyiapws95/

License
-------
pyiapws95 is is released under
`MIT Licence <https://choosealicense.com/licenses/mit>`__.
