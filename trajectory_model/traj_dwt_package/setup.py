#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

# Get version from __init__.py
with open('src/traj_dwt/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__version__'):
            version = line.split('=')[1].strip().strip('"\'')
            break

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='traj_dwt',
    version=version,
    description='Trajectory Dynamic Time Warping Package for gene expression analysis',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Author Name',
    author_email='author@example.com',
    url='https://github.com/yourusername/traj_dwt',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.19.0',
        'scipy>=1.5.0',
        'pandas>=1.0.0',
        'matplotlib>=3.3.0',
        'scikit-learn>=0.23.0',
        'fastdtw>=0.3.0',
        'joblib>=0.16.0',
    ],
    extras_require={
        'dev': [
            'pytest>=6.0.0',
            'pytest-cov>=2.10.0',
            'flake8>=3.8.0',
            'black>=20.8b1',
        ],
        'anndata': [
            'anndata>=0.7.6',
            'scanpy>=1.6.0',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
    keywords='trajectory, gene expression, dynamic time warping, bioinformatics, single-cell',
    project_urls={
        'Documentation': 'https://github.com/yourusername/traj_dwt',
        'Source': 'https://github.com/yourusername/traj_dwt',
        'Tracker': 'https://github.com/yourusername/traj_dwt/issues',
    },
) 