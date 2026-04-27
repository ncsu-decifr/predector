from setuptools import find_packages
from setuptools import setup

setup_requires = [
    'numpy',
]

install_requires = [
    'numpy',
    'scipy',
    'Theano==1.0.1',
    'Lasagne==0.2.dev1'
    ]

setup(
    name='DeepLoc',
    version='0.1',
    author='Jose Juan Almagro Armenteros',
    author_email='jjaa@bioinformatics.dtu.dk',
    packages=find_packages(),
    scripts=['bin/deeploc'],
    setup_requires=setup_requires,
    install_requires=install_requires,
    dependency_links=['https://github.com/Lasagne/Lasagne/tarball/master#egg=Lasagne-0.2.dev1'],
    package_data={'DeepLoc': ['parameters/*.npz']},
    description='Prediction of protein subcellular localization for eukaryotic proteins.',
    long_description=open('README.txt').read(),
)
