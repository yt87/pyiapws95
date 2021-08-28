'''
setup.py for pyiapws95
'''

from setuptools import setup
import versioneer

CLASSIFIERS = [
    "Development Status :: 3 - Alpha",
    "License :: Public Domain",
    "Operating System :: OS Independent",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Topic :: Scientific/Engineering",
]

metadata = dict(
    name="pyiapws95",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Python package to compute thermodynamic properties of water.",
    author='George Trojan',
    author_email='george.trojan@gmail.com',
    license='BSD-0-Clause',
    url="https://github.com/yt87/pyiapws95",
    install_requires=["numpy>=1.18"],
    extras_require={"pint": "pint>=0.10", "numba": "numba>=0.48"},
    classifiers=CLASSIFIERS,
    test_suite="tests",
    tests_require=["pytest"],
    zip_safe="True",
    )

setup(**metadata)
