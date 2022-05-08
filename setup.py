import os
from setuptools import setup


def version():
    setup_dir = os.path.dirname(os.path.realpath(__file__))
    version_file = open(os.path.join(setup_dir, 'MarkerMAG', 'VERSION'))
    return version_file.readline().strip()


__long_description__ = '''
MarkerMAG: linking MAGs with 16S rRNA marker genes

Weizhi Song, Torsten Thomas

Center for Marine Science & Innovation (CMSI)

University of New South Wales, Sydney, Australia
'''


setup(name="MarkerMAG",
      description="linking MAGs with 16S rRNA marker genes",
      keywords="Bioinformatics Metagenomics MAG 16S",
      author="Weizhi Song, Torsten Thomas",
      author_email="songwz03@gmail.com",
      version=version(),
      long_description=__long_description__,
      license="GPL3+",
      #url="https://github.com/songweizhi/MarkerMAG",
      url="https://pypi.org/project/MarkerMAG/",
      # classifiers=[
      #       'Development Status :: 5 - Production/Stable',
      #       'Intended Audience :: Science/Research',
      #       'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      #       'Natural Language :: English',
      #       'Programming Language :: Python :: 3.10',
      #       'Topic :: Scientific/Engineering :: Bio-Informatics'],
      scripts=['bin/MarkerMAG'],
      packages=['MarkerMAG', 'MarkerMAG.vxtractor'],
      package_data={'MarkerMAG': ['VERSION', 'SILVA_16S_order.fasta'],
                    'MarkerMAG.vxtractor': ['gpl.txt',
                                            'vxtractor.pl',
                                            'readmefungiLSU.txt',
                                            'HMMs/SSU_archaea/*',
                                            'HMMs/SSU_bacteria/*',
                                            'README.txt']},
      include_package_data=True,
      install_requires=['biopython', 'pandas', 'plotly', 'seaborn', 'matplotlib', 'numpy', 'setuptools'])
      # zip_safe=False
