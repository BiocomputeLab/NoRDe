#!/usr/bin/env python3
"""
Setup script for NoRDe (Nonrepetitive RNA Designer)
"""

from setuptools import setup, find_packages
import os

# Read the README file
def read_readme():
    with open("README", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="norde",
    version="1.1.0",
    author="Konstantinos Kariotis",
    author_email="kariotis.konst@gmail.com",
    description="A modular toolkit for RNA variant design and analysis",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="https://github.com/kkariotis/rna-variant-tool",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.0",
            "black>=21.0",
            "flake8>=3.8",
            "mypy>=0.800",
        ],
        "docs": [
            "sphinx>=4.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "norde=tool.cli:run_cli",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    keywords="rna, bioinformatics, variant-design, structure-prediction, conservation-analysis",
    project_urls={
        "Bug Reports": "https://github.com/kkariotis/rna-variant-tool/issues",
        "Source": "https://github.com/kkariotis/rna-variant-tool",
        "Documentation": "https://github.com/kkariotis/rna-variant-tool#readme",
    },
) 