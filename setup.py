try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

version = '0.1.3.dev0'
setup(
    name='emstimtools',
    version=version,
    author='Julius Zimmermann',
    author_email='julius.zimmermann@uni-rostock.de',
    packages=['emstimtools', 'emstimtools.utils', 'emstimtools.fenics', 'emstimtools.dataanalysis'],
    description='EMStimTools is a package that provides a SALOME-GMSH-FEniCS workflow to solve problems related to electromagnetic stimulation of for instance biological tissue or cell cultures.',
    install_requires=['pyyaml==5.3'],
    long_description=open('README.rst').read()
)
