try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'pyrxn',
    'author': 'Justin Vrana',
    'url': '',
    'download_url': '',
    'author_email': 'justin.vrana@gmail.com',
    'version': '0.0.2',
    'install_requires': ['pandas', 'numpy', 'scipy'],
    'packages': ['pyrxn'],
    'scripts': [],
    'name': 'pyrxn',
    'license': 'Copyright University of Washington'
}

setup(**config)
