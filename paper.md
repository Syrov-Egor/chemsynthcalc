---
title: 'chemsynthcalc: A Software For Chemical Synthesis Calculations And Reaction Balancing'
tags:
  - Python
  - chemistry
  - chemical reactions
  - linear algebra
authors:
  - name: Egor V. Syrov
    orcid: 0000-0002-6394-2367
    affiliation: 1
affiliations:
 - name: N. I. Lobachevsky State University of Nizhny Novgorod, Nizhny Novgorod, Russia
   index: 1
   ror: 01bb1zm18
date: 19 January 2026
bibliography: paper.bib
---

# Summary

The problem of determining the amounts of reactants needed for a chemical 
synthesis (and the related problem of balancing reaction equations) 
is still relevant. In this article, we present `chemsynthcalc`, 
a Python package and a GUI interface for automatic synthesis calculations 
and reaction balancing. The main advantages of this package over its 
competitors are operations with non-integer coefficients and atomic amounts, 
and a fast, rich, and well-documented API. Four linear algebra methods are implemented for 
performing operations on a chemical reaction matrix. Three of them are 
adopted from the literature and one of them (the combinatorial or exhaustive search method) is 
original. The stability and performance of the calculations were 
measured against a large dataset of inorganic synthesis reactions.

# Statement of need

The solid-state synthesis of inorganic compounds is a vital field in 
material science for the production of a variety of materials. The synthesis
of photocatalysts, superconductors, ferroelectrics, phosphors, and other materials
by calcination of solids, hydrothermal or sol-gel methods requires 
calculations of the precursor masses before weighing, grinding, and 
heating.

`chemsynthcalc` was created to help scientists calculate inorganic synthesis. 
There are three main aspects to this problem. First, a single student, 
scientist, or postdoc, working in their lab trying to synthesize a new material. 
The sample size of potential reactions is quite small in this case, 
and the scientist can use the traditional pen-and-paper method. 
However, automation of this routine task offers clear advantages.
Second, while we are not quite there yet, one can imagine a 
robotic inorganic synthesis station, such as those developed by Dr. Cronin and colleagues [@Sans2015].
In such case, a large number of reactions must be calculated quickly and precisely,
and we really should not hard-code all the masses beforehand. 
Finally, there are datasets of text-mined inorganic reactions [@Kononova2019] that will 
surely expand and grow in size. With a sample size of thousands, 
reaction balancing software must be extremely fast, robust, and flexible enough to balance 
as many reactions as possible. In this case, a free, open-source solution that can be 
embedded in data processing pipelines is essential.

`chemsynthcalc` addresses all three of those cases. 
It is simple enough to use by a single scientist 
familiar with using Python packages, yet fast and robust enough to 
precisely calculate hundreds and thousands of reactions. There are 
already a large number of reaction balancing software, 
why do one need `chemsynthcalc` and why is it better than competitors? 
Here are some advantages of this package:

* It is completely free and open-source (under MIT license).

* It provides a rich and simple API for its functions.

* It can balance a huge variety of reactions.

* It can deal with formulas with float atom counts
  (such as RbLa~0.99~Eu~0.01~Nb~2~O~7~), both for molar mass 
  calculations and reaction balancing.

* It supports a large number of nested parentheses
  in formulas (up to recursion depth limit).

* It supports addition notation (like CuSO~4~\*5H~2~O).

* It is fast due to NumPy matrix operations.

* Finally, it is one of the few programs that 
  can directly output precursor masses from a reaction 
  string in three lines of code.

# Software design

The design philosophy of `chemsynthcalc` is based on creating a 
friendly, yet deep and powerful API for users of any level 
of programming proficiency. Most of the other reaction balancing
and mass calculation software are closed source and are offered 
as a service (such as WebQC [@webqc_org_balance]). While they can be 
suitable for infrequent use, they are still lacking all the 
aforementioned benefits of FOSS. The only mature open-source Python 
package with similar capabilities is `ChemPy` [@Dahlgren2018]. This 
package is using `sympy.linsolve()` function to symbolically 
solve a system of linear equations. While the symbolic method is 
suitable for simple reactions (and can even perform some operations 
that `chemsynthcalc` cannot, such as balancing reactions with ionic 
and abstract formulas), it is lacking precision and rigor required for complex 
reactions with non-integer atomic amounts. The other approach is a 
strictly numeric linear algebra reaction balancing with NumPy, which was 
chosen for `chemsynthcalc`. Therefore, contributing a complete feature 
overhaul to the `ChemPy` package would be inappropriate. Moreover, 
`ChempPy` is a large package with main focus on equlibria and chemical 
kinetics, while we prefer more modular approach with (almost) 
single-purpose software with minimal external dependencies.

# Research impact statement

While the `chemsynthcalc` package is not widely adopted yet, it is extensively 
used by students of different levels in the Lobachevsky State University of 
Nizhny Novgorod. With user-friendly web, desktop, and Android applications it is helpful 
for day-to-day use in any inorganic chemistry lab. It is can also be used 
by chemistry teachers and students to study chemical reactions at all 
levels of education. Moreover, in the future it could be used to fully 
automate one of the steps of inorganic synthesis. It is also ready for use in scanning 
huge datasets of reactions for balancing and mass calculation purposes.

`chemsynthcalc` has approximately 24,000 all-time downloads via PyPI. The comprehensive unit 
test suite (96% test coverage), along with detailed automatic documentation, 
modular structure and CI script, makes it easy to contribute new features to the project 
(like single line solid solution reactions support, isotopes, etc.). 
`chemsynthcalc` is also benchmarked against a dataset of real-world data-mined 
solid-state reactions (see Testing section).

# Algorithm

The calculation can be described algorithmically as follows:

1. Selection of synthesis *target* and *target mass*.
  (which substance will be synthesized and in what quantity).

2. Chemical equation construction (while taking available 
  precursors into account).

3. Chemical equation balancing.

4. Calculation of molar masses of compounds.

5. Calculation of the amount of substance of the target compound
  *n* (mol) as $n=m/M$, where *m* is target mass 
  (in grams), and *M* is the target molar mass (in g/mol).

6. Broadcasting this amount of substance to other compounds
  considering their stoichiometry and calculation of their masses as
  $m=cnM$ (where *c* is a ratio of stoichiometric coefficients
  of a given substance and target substance).

Two most challenging steps are reaction balancing 
and calculation of molar mass. There is a set of algebraic or mathematical balancing methods.
They are generally based on creation of a system of linear equations 
for every atom in the chemical equation and solving this system using 
different linear algebra techniques. These 
methods are based on the concept of chemical composition matrix or 
chemical reaction matrix [@Blakley1982]. For example, the chemical reaction matrix
for a simple textbook reaction
$$
KMnO_4 + HCl = MnCl_2 + Cl_2 + H_2O + KCl
$$
will look like:
$$
\begin{matrix}
    K \\
    Mn \\
    O \\
    H \\
    Cl \\
\end{matrix}
\begin{bmatrix}
    1 & 0 & 0 & 0 & 0 & 1 \\
    1 & 0 & 1 & 0 & 0 & 0 \\
    4 & 0 & 0 & 0 & 1 & 0 \\
    0 & 1 & 0 & 0 & 2 & 0 \\
    0 & 1 & 2 & 2 & 0 & 1 \\
\end{bmatrix}
$$
where rows represent the atom type, and columns represent the compound.
Now, this system of linear equations can be solved using linear algebra techniques.
In `chemsynthcalc`, there are four such methods: the `inverse` method based on Thorne paper
[@Thorne2011]; the `general pseudoinverse` and `partial pseudoinverse` methods based on 
Risteski papers [@Risteski2008; @Risteski2009; @Risteski2013]. Finally, the `combinatorial` 
method is used to solve Diophantine matrix equation $Ax=By$ by brute force. 
Despite its limitations, the combinatorial method can achieve some unexpected and 
interesting results, for example, finding a smaller set of coefficients for a well-known 
reaction. Detailed descriptions of these algorithms are presented in the documentation.

There is also a user-friendly, full-featured, cross-platform GUI version of this package 
built with a Go backend and web frontend, available for Windows, Linux, macOS, 
Android, and web [@chemsynthcalc-GUI].

# Testing

Beside standard unit tests, it is important to test and benchmark such a package against 
a wide variety of real-world examples of reactions. A dataset of text-mined solid-state reactions 
[@Kononova2019] was used to test its capabilities. The reaction 
dataset used was a part of the `solid-state_dataset_20200713` file. 
The original dataset contains more than 30,000 entries. We filtered this list of reactions, 
retaining only valid deduplicated reactions. These do not include 
reactions that contains non-stoichiometry symbols (like $\Delta$), 
reactions with unknown amounts of substance (like $x$ or $y$), or 
reaction with different sets of atoms on right and left sides of 
the equation.

Thus, the list of 9,181 valid reactions was formed. 
Then automatic balancing and calculations of the masses were performed 
for each reaction on this list. The results (benchmarked on Ubuntu 24.04, 
AMD Ryzen 7 5700x, and 64 GB DDR4 RAM with Python 3.13.9, NumPy 2.4.0 
and output to a .txt file) are 0.057 ms per formula and 0.49 ms per reaction. 
There are 69 million reactions in the Reaxys database [@ReaxysDatabase]; 
therefore all the reactions in the largest database can be balanced within
9.5 hours on a single consumer-grade PC.

# AI usage disclosure

No generative AI tools were used in the development of this software or the writing
of this manuscript.

# References
