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

VERSION = '1.0.1'

config = {
    'description': 'A lattice-gas Monte-Carlo simulation tool',
    'long_description': read('README.md'),
    'author': 'Benjamin J. Morgan',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/lattice_mc',
    'download_url': "https://github.com/bjmorgan/lattice_mc/archive/%s.tar.gz" % (VERSION),
    'author_email': 'b.j.morgan@bath.ac.uk',
    'version': VERSION,
    'install_requires': [ 'numpy', 
                          'matplotlib', 
                          'pandas', 
                          'scipy', 
                          'coverage==4.3.4',
                          'codeclimate-test-reporter' ],
    'license': 'MIT',
    'packages': ['lattice_mc'],
    'scripts': [],
    'name': 'lattice_mc'
}

setup(**config)
