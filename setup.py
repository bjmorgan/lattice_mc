import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

config = {
    'description': 'A lattice-gas Monte-Carlo simulation tool',
    'long_description': read('README'),
    'author': 'Benjamin J. Morgan',
    'url': 'https://github.com/bjmorgan/lattice_mc',
    'download_url': 'https://github.com/bjmorgan/lattice_mc/archive/master.zip',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'version': '0.1',
    'install_requires': ['numpy'],
    'license': 'MIT',
    'packages': ['lattice_mc', 'tests'],
    'scripts': [],
    'name': 'lattice_mc'
}

setup(**config)
