import os
from setuptools import setup


def version():

    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'MarkerMAG', 'VERSION'))

    return version_file.readline().strip()


__long_description__ = '''

MarkerMAG: link MAGs with marker genes

Weizhi Song (songwz03@gmail.com)

Center for Marine Science & Innovation (CMSI)

University of New South Wales, Sydney, Australia

'''

setup(name="MarkerMAG",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      author="Weizhi Song, Torsten Thomas",
      author_email="songwz03@gmail.com",
      keywords="Bioinformatics Metagenomics Marker_Gene",
      description="link MAGs with marker genes",
      url="https://github.com/songweizhi/MarkerMAG",
      packages=['MarkerMAG'],
      package_data={'': ['*.r', '*.R', '*.py', 'VERSION']},
      include_package_data=True,
      install_requires=['biopython', 'pandas'],
      scripts=['bin/MarkerMAG'])
