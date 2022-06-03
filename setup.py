from pathlib import Path
from setuptools import setup

this_directory = Path(__file__).parent;
long_description = (this_directory / "README.md").read_text();

VERSION = '0.0.1';
DESCRIPTION = 'Herramientas para análisis de secuencias en resultados de ChIP-seq.';
PACKAGE_NAME = 'analisischip';
AUTHOR = 'Emilio Kolomenski';
EMAIL = 'ekolomenski@gmail.com';
GITHUB_URL = 'https://github.com/EmilioKolo/analisissecuenciaschip';

setup(
    name = PACKAGE_NAME,
    packages = [PACKAGE_NAME],
    version = VERSION,
    license='MIT',
    description = DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    author = AUTHOR,
    author_email = EMAIL,
    url = GITHUB_URL,
    keywords = [''],
    install_requires=[ 
        'numpy', 
        'matplotlib', 
        'pyensembl'
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
)