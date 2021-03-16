from skbuild import setup
from setuptools import find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="spyprot",
    version='0.3.5',
    author="INTERDISCIPLINARY LABORATORY of BIOLOGICAL SYSTEMS MODELLING, University of Warsaw, Warsaw, Poland",
    author_email="bmjastrzebski@gmail.com, p.rubach@cent.uw.edu.pl",
    description="This package provides a set of tools for accessing protein databases and manipulating PDB/CIF files.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ilbsm/spyProt",
    packages=find_packages(),
    install_requires=[
        'tqdm>=4.31.1', 'lxml>=4.5.0', 'requests>=2.20.0', 'biopython==1.76', 'psutil>5.6.6', 'subprocess32>=3.5.0', 'wget>=3.2', 'mysolr>=0.8.3'
    ],
    python_requires='>=2.7.0',
    classifiers=[
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"
    ],
)
