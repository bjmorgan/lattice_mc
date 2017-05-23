import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst')
except ImportError:
    long_description = open('README.md').read()

VERSION = '1.0.2'

config = {
    'description': 'A lattice-gas Monte-Carlo simulation tool',
    'long_description': long_description,
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
