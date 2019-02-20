import os
import re
from distutils.core import setup


tests_require = [
    'pytest',
    'pytest-cov',
    'pytest-benchmark'
]

install_requires = [
    'pandas',
    'numpy',
    'scipy',
    'seaborn',
]

NAME = 'pyrxn'
MAIN_PKG = 'pyrxn'

def parse_version_file():
    here = os.path.abspath(os.path.dirname(__file__))
    ver_dict = {}
    with open(os.path.join(here, MAIN_PKG, '__version__.py'), 'r') as f:
        for line in f.readlines():
            m = re.match('__(\w+)__\s*=\s*(.+)', line)
            if m:
                ver_dict[m.group(1)] = m.group(2)
    return ver_dict


ver = parse_version_file()


setup(
    title=ver['title'],
    name=NAME,
    version=ver['version'],
    packages=[MAIN_PKG],
    url=ver['url'],
    license='MIT',
    author=ver['author'],
    author_email='justin.vrana@gmail.com',
    description='',
    install_requires=install_requires,
    tests_require=tests_require,
)
