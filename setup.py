import pkg_resources, os, pkgutil
from setuptools import setup, find_packages
from version import __version__
import glob

def dependencies():
    file_ = pkg_resources.resource_filename(__name__, os.path.join('requirements', 'default.txt'))
    with open(file_, 'r') as f:
        return f.read().splitlines()

setup(
    name                 = 'amber',
    version              = __version__,
    description          = 'AMBER: Assessment of Metagenome BinnERs',
    author               = 'CAMI',
    author_email         = 'contact@cami-challenge.org',
    url                  = 'http://cami-challenge.org',
    scripts              = glob.glob('*.py'),
    install_requires     = dependencies(),
    packages             = find_packages(),
    classifiers = [
        'Natural Language :: English',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.5',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX'
    ]
)
