#!/usr/bin/env python
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="zonkey",
    version="0.1.0",
    author="Eliot Boulanger",
    author_email="eliot.boulanger@kuleuven.be",
    description="An open source QM/MM package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/eliotblg/zonkey.git",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3",
    ),
)
