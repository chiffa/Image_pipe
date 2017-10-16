"""
Generic Setup configured for the project
"""
from setuptools import setup, find_packages  # Always prefer setuptools over distutils
from codecs import open  # To use a consistent encoding
from os import path
import os

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

requirements = [
    'pillow',
    'numpy',
    'scipy',
    'scikit-image',
    'matplotlib',
    ]

setup(
    name='ImagePipe',
    version='0.2.0_rc1',
    description='Biological image analysis pipelines compiler',
    long_description=long_description,
    url='https://github.com/chiffa/Image_pipe',
    author='Andrei Kucharavy',
    author_email='andrei.shield15@gmail.com',
    license='BSD',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research, Developers',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Image Recognition',
        'Topic :: Scientific/Engineering :: Information analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Microsoft :: Windows :: Windows 10',
        'Operating System :: Microsoft :: Windows :: Windows 7',
        'Operating System :: Microsoft :: Windows :: Windows 8',
        'License :: OSI Approved :: BSD 3-clause license',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: CPython',
    ],
    keywords='image analysis, 4D microscopy, image analysis pipeline',
    packages=['imagepipe']+find_packages(),
    install_requires=requirements,
)