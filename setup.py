#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['numpy',
                'xarray',
                'netcdf4',
                'pandas',
                'scikit-image',
                'scipy',
                'tqdm',
                'shapely',
                'scikit-learn',
                'matplotlib',
                'cartopy']

test_requirements = ['pytest>=3', ]

setup(
    author="Severin Kaderli",
    author_email='severin.kaderli@unibe.ch',
    python_requires='>=3.8',
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
	    'Intended Audience :: Science/Research',
	    'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],
    description="Detect, classify, and track Rossby Wave Breaking (RWB).",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/x-rst',
    include_package_data=True,
    keywords='wavebreaking',
    name='wavebreaking',
    packages=find_packages(include=['wavebreaking', 'wavebreaking.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/skaderli/wavebreaking',
    version='0.3.0',
    zip_safe=False,
)
