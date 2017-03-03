import os
import sys
from setuptools import setup, find_packages
from setuptools.extension import Extension

try:
    import numpy as np
except ImportError:
    raise ImportError("numpy must be installed prior to installing biom")

setup(name='denovo_picrust',
      version='0.1',
      description='Pipeline to run PICRUSt on any OTU picking method',
      long_description='This pipeline uses ancestral state reconstruction\
        to allow you to run PICRUSt on any "OTUs"',
      classifiers=[
        'Development Status :: Alpha',
        'License :: GNU GPL V3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Bioinformatics :: Microbiome :: Metagenomics',
      ],
      keywords='Bioinformatics Microbiome Metagenomics',
      url='http://github.com/vmaffei/dada2_to_picrust',
      author='Gene Blanchard, Vince Maffei',
      maintainer='Gene Blanchard',
      license='GNU GPL V3',
      scripts=['bin/denovo_picrust.py','bin/gap_trimmer.py'],
      include_package_data=True,
      include_dirs=[np.get_include()],
      package_data={'denovo_picrust': ['data/*']},
      install_requires=[
          'numpy==1.7.1',
          'pynast==1.2.2',
      ],
      zip_safe=True)
