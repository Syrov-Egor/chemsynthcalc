from setuptools import setup

setup(
   name='chemsynthcalc',
   version="1.0.4",
   description="Python package for calculating the masses of substances required for chemical synthesis directly from the reaction string. It includes solutions for all intermidiate steps, including chemical formula parsing, molar mass calculation and reaction balancing with different matrix methods.",
   author="Egor Syrov",
   author_email="syrov@unn.ru",
   packages=['chemsynthcalc'],  #same as name
   install_requires=["scipy>=1.6"], #external packages as dependencies
)