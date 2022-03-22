#! /usr/bin/env python

import pathlib
from setuptools import setup, find_packages

__here__ = pathlib.Path(__file__).parent

extra_files = {"phasematching_calc": ["VERSION"]}

with __here__ / "phasematching_calc" / "VERSION" as version_file:
    version = version_file.read_text().strip()

setup(
    name="phasematching_calc",
    packages=find_packages(),
    package_data=extra_files,
    install_requires=["numpy", "WrightTools"],

    version=version,
    description="A simulation package for simulating phasematching effects in discrete layer isotropic samples.",
    author="phasematching_calc Developers",
    license="MIT",
    url="https://github.com/wright-group/phasematching_calc",
    keywords="spectroscopy science multidimensional simulation",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
    ],
)
