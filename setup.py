#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

requirements = [
    "numpy",
    "hanselx",
    "pysam",
    "PyVCF",
]

test_requirements = [

]

setuptools.setup(
    name="gretel",
    version="0.0.8",
    url="https://github.com/samstudio8/gretel",

    description="An algorithm for recovering potential haplotypes from metagenomes",
    long_description="",

    author="Sam Nicholls",
    author_email="sam@samnicholls.net",

    maintainer="Sam Nicholls",
    maintainer_email="sam@samnicholls.net",

    packages=setuptools.find_packages(),
    include_package_data=True,

    install_requires=requirements,

    entry_points = {
        "console_scripts": [
            "gretel=gretel.cmd:main",
        ]
    },

    classifiers = [
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],

    test_suite="tests",
    tests_require=test_requirements
)
