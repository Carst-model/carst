#!/usr/bin/env python

from distutils.core import setup
from glob import glob

import versioneer

CMD_CLASS = versioneer.get_cmdclass()

setup(
    name='carst',
    cmdclass=CMD_CLASS,
    version=versioneer.get_version(),
    description='Finite element stratigraphic model',
    author='Jon Hill',
    author_email='jon.hill@york.ac.uk',
    url='https://github.com/Carst-model/carst',
    packages=[
        'carst',
        'test',
        'examples',
        'carst_config',
    ],
    scripts=glob('scripts/*'),
)
