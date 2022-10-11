.. chemsynthcalc documentation master file, created by
   sphinx-quickstart on Sun Sep 11 16:13:47 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

chemsynthcalc
=============
**ChemSynthCalc** stands for **Chem**\ ical **Synth**\ esis **Calc**\ ulator - 
a Python 3 package for calculating the masses of substances required for 
chemical synthesis directly from the reaction string.
It includes solutions for all intermidiate steps:

* Chemical formula strings parsing
* Molar mass calculation
* Reaction string parsing and decomposition
* Reaction matrix creation
* Manual reaction balance
* Automatic reaction balance
* Masses calculations
* Output to terminal, as plain txt or JSON (object or file)
* ... and more

Background
==========
ChemSynthCalc was created to address the issue of inorganic synthesis calculations.
There are three main aspects to this problem. 

First of all, a case of single student, scientist, postdoc, etc. sitting in their lab 
and trying to synthesize a new material. The sample size of potential reactions is 
quite small in this case, and our fellow scientist can use an the old-fashioned way 
of pen and paper. But why one should not try to automate this boring task?

Secondly, while we are not quite there yet, one can imagine a robotic inorganic synthesis
station, like `Dr Cronin <https://pubs.rsc.org/en/content/articlelanding/2012/LC/c2lc40761b>`_
devices. In this case, we need to calculate a large number of reactions fast and precise, and
we really should not hardcore all the masses beforehand.

Finally, there are datasets of `text-mined inorganic reaction <https://www.nature.com/articles/s41597-019-0224-1>`_.
They surely will expand and grow in size. With sample size of ten of thousands we need our reaction balancing
software to be extremely fast, robust yet flexible enough to balance as many reactions as it can. In this case,
we surely need a free, open-source solution that can be embedded in the data processing system.

ChemSynthCalc addresses all three of those cases. It is simple enough to use for a single
scientist who are familiar with using Python packages, and fast and robust enough to precisely
calculate hundreds and thousands reactions.

Contents
========
.. toctree::
   :maxdepth: 3

   installation
   usage
   chemical_formula
   chemical_reaction
   API
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`