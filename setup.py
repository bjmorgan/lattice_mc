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
    'long_description': read('README.md'),
    'author': 'Benjamin J. Morgan',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/lattice-mc',
    'download_url': 'https://github.com/bjmorgan/lattice-mc/tarball/0.9.1',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'version': '1.0.0',
    'install_requires': ['numpy', 'matplotlib', 'pandas'],
    'license': 'MIT',
    'packages': ['lattice-mc'],
    'scripts': [],
    'name': 'lattice-mc'
}

setup(**config)
