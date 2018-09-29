import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

long_description = open('README.md').read()

from lattice_mc import __version__
VERSION = __version__

config = {
    'description': 'A lattice-gas Monte-Carlo simulation tool',
    'long_description': long_description,
    'long_description_content_type': 'text/markdown',
    'author': 'Benjamin J. Morgan',
    'author_email': 'b.j.morgan@bath.ac.uk',
    'url': 'https://github.com/bjmorgan/lattice_mc',
    'download_url': "https://github.com/bjmorgan/lattice_mc/archive/%s.tar.gz" % (VERSION),
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
