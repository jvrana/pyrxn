try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
	from distutils.command.build_py import build_py_2to3 \
		as build_py
except ImportError:
	from distutils.command.build_py import build_py

config = {
    'description': 'pyrxn',
    'author': 'Justin Vrana',
    'url': '',
    'download_url': '',
    'author_email': 'justin.vrana@gmail.com',
    'version': '0.0.1',
    #'install_requires': ['pandas', 'numpy', 'scipy'],
    'packages': ['pyrxn'],
    'scripts': [],
    'name': 'pyrxn',
    'license': 'Copyright University of Washington'
}

setup(**config)
