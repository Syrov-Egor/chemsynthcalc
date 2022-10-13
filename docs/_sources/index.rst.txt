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
* Calculation of molar mass
* Reaction string parsing and decomposition
* Creation of reaction matrix
* Manual reaction balance
* Automatic reaction balance
* Masses calculations
* Output to terminal, as plain txt or JSON (object or file)
* ... and more

Main project page
=================
https://github.com/Syrov-Egor/chemsynthcalc

Background
==========
ChemSynthCalc was created to address the issue of inorganic synthesis calculations.
There are three main aspects to this problem. 

First of all, a single student, scientist, postdoc, etc., sitting in their lab 
trying to synthesize a new material. The sample size of potential reactions is quite small in this case, 
and our fellow scientist can use the old-fashioned pen and paper method. 
But why should one not try to automate this boring task?

Secondly, while we are not quite there yet, one can imagine a robotic inorganic synthesis
station, like `Dr. Cronin <https://pubs.rsc.org/en/content/articlelanding/2012/LC/c2lc40761b>`_
devices. In this case, we need to calculate a large number of reactions fast and precisely, 
and we really should not hard-core all the masses beforehand.

Finally, there are datasets of `text-mined inorganic reactions <https://www.nature.com/articles/s41597-019-0224-1>`_.
They will surely expand and grow in size. With a sample size of tens of thousands, we need our reaction balancing software 
to be extremely fast, robust, and flexible enough to balance as many reactions as it can. 
In this case, we surely need a free open-source solution that can be embedded in the data processing system.

ChemSynthCalc addresses all three of those cases. It is simple enough to use for a single scientist who 
is familiar with using Python packages and fast and robust enough to precisely 
calculate hundreds and thousands of reactions.

Competitive analysis
======================
There are already a large number of reaction balancing software, why do you need
chemsynthcalc and why is it better than competitors?

* chemsynthcalc is completely free and open-source (under MIT license)
* chemsynthcalc provides a rich and simple API for its functions
* chemsynthcalc can balance a huge variety of reactions
* chemsynthcalc can deal with formulas with float atom count (like RbLa0.99Eu0.01Nb2O7), both for molar mass and reaction balancing
* chemsynthcalc supports an infinite amount of nested parentheses in formulas
* chemsynthcalc supports adduct notaion (like CuSO4*5H2O)
* chemsynthcalc is fast thanks to NumPy matrix operations
* finally, chemsynthcalc is one of the few programs that can directly output precursor masses from a reaction string (in three lines of code!)

.. list-table:: chemsynthcalc competitive analysis
   :widths: 25 25 25 25
   :header-rows: 1

   * - Software
     - Mass calculation
     - Float coefficients and atom amounts
     - API
   * - `chemsynthcalc <https://github.com/Syrov-Egor/chemsynthcalc>`_
     - ✓
     - ✓
     - ✓
   * - `WebQC <https://www.webqc.org/balance.php>`_
     - ✓
     - x
     - x
   * - `Equation Balancer <https://equationbalancer.com/>`_
     - x
     - x
     - x
   * - `ChemicalAid <https://en.intl.chemicalaid.com/tools/equationbalancer.php>`_
     - x
     - ✓
     - x
   * - `CHEMIX School <https://www.chemix-chemistry-software.com/chemistry-software.html>`_
     - x
     - ✓
     - x
   * - `ChemPy <https://github.com/bjodah/chempy>`_
     - x
     - x
     - ✓

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
   contacts


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`