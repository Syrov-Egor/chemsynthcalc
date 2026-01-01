---
title: 'chemsynthcalc: a program for balancing and calculating chemical reactions'
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
date: 01 January 2026
bibliography: paper.bib
---

# Summary

The problem of finding the amounts of reactants needed for a chemical 
synthesis (and the related problem of balancing reaction equations) 
is still relevant. In this article, we present `chemsynthcalc` - 
a Python package and GUI interface for automatic synthesis calculations 
and reaction balancing. The main advantages of this package over its 
competitors are operations with non-integer coefficients and atomic amounts, 
and a fast, rich, and well-documented API. Four linear algebra methods are implemented for 
performing operations on a chemical reaction matrix. Three of them are 
adopted from the literature, and one of them (combinatorial or exhaustive search method) is 
original. The stability and performance of the calculations were 
measured against a large dataset of inorganic synthesis reactions.

# Statement of need

A solid-state synthesis of inorganic compounds is vital field in 
material science for production of variety of materials. The synthesis
of photocatalysts, superconductors, ferroelectrics, phosphors and other materials
by calcination of solids, hydrothermal or sol-gel methods requires
calculations of the precursor masses before weighing, grinding and 
heating.

`chemsynthcalc` was created to help scientists calculate inorganic synthesis. 
There are three main aspects to this problem. First of all, a single student, 
scientist, postdoc, etc., working in their lab trying to synthesize a new material. 
The sample size of potential reactions is quite small in this case, 
and the scientist can use the old-fashioned pen-and-paper method. 
But why should one not try to automate this boring task? 
Secondly, while we are not quite there yet, one can imagine a 
robotic inorganic synthesis station, like Dr. Cronin [@Sans2015] devices. 
In this case, we need to calculate a large number of reactions fast and precisely,
and we really should not hard-code all the masses beforehand. 
Finally, there are datasets of text-mined inorganic reactions [@Kononova2019]. 
They will surely expand and grow in size. With a sample size of thousands, 
we need our reaction balancing software to be extremely fast, robust, and 
flexible enough to balance as many reactions as it can. In this case, 
we need a free open-source solution that can be embedded in 
the data processing pipeline.

`chemsynthcalc` addresses all three of those cases. 
It is simple enough to use for a single scientist who is 
familiar with using Python packages, and fast and robust enough to 
precisely calculate hundreds and thousands of reactions. There are 
already a large number of reaction balancing software, 
why do one may need `chemsynthcalc` and why is it better than competitors? 
Here are some advantages of this package:

* It is completely free and open-source (under MIT license).

* It provides a rich and simple API for its functions.

* It can balance a huge variety of reactions.

* It can deal with formulas with float atom count 
  (like RbLa~0.99~Eu~0.01~Nb~2~O~7~), both for molar mass and 
  reaction balancing.

* It supports an infinite number of nested parentheses
  in formulas.

* It supports addition notation (like CuSO~4~\*5H~2~O).

* It is fast thanks to NumPy matrix operations.

* Finally, it is one of the few programs that 
  can directly output precursor masses from a reaction 
  string (in three lines of code!).

# Algorithm

The calculation can be algorithmized as following:

1. Selection of synthesis *target* and *target mass*.
  (which substance will be synthesized and in what quantity).

2. Chemical equation construction (while taking available 
  precursors into account).

3. Chemical equation balancing.

4. Calculation of molar masses of compounds.

5. Calculation of the amount of substance of the target compound
  *n* (mol) as $(n=m/M)$, where *m* is target mass 
  (in grams), and *M* is the target molar mass (in g/mol).

6. Broadcasting this amount of substance to other compounds
  considering their stoichiometry and calculation of their masses as
  $(m=cnM)$ (where *c* is a ratio of stoichiometric coefficients
  of some substance and target substance).

Two most challenging steps here are, of course, the reaction balancing 
and the calculation of the molar mass. There is a set of algebraic or mathematical balancing methods.
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
Where rows represent the atom type, and columns represent the compound.
Now, this system of linear equations can be solved using linear algebra techniques.
In `chemsynthcalc`, there are four such methods: the `inverse` method based on Thorne paper
[@Thorne2011]; the `general pseudoinverse` and `partial pseudoinverse` methods based on 
Risteski papers [@Risteski2008; @Risteski2009; @Risteski2013]. Finally, the `combinatorial` 
method is used to solve Diophantine matrix equation $Ax=By$ by brute force. 
Despite its limitations, the combinatorial method can achieve some unexpected and 
interesting results, for example, find a smaller set of coefficients for a well known 
reaction. Detailed descriptions of this algorithms are presented in the docs.

There is also a user-friendly full-functional crossplatform GUI version of this package, 
available for Windows, Linux, MacOS, Android and web [@chemsynthcalc-GUI].

# Testing
Beside standard unit tests, it is important to test and benchmark such package against 
huge variety of real-life examples of reactions. A dataset of text-mined solid-state reactions 
[@Kononova2019] was used to test its capabilities. The reaction 
dataset used was a part of the `solid-state_dataset_20200713` file. 
The original dataset contains more than 30,000 entries on inorganic 
synthesis, including reactions. We filter this list of reactions, 
leaving only valid deduplicated reactions. These do not include 
a reactions that contains non-stoichiometry symbols (like $\delta$), 
reactions with unknown amounts of substance (like $x$ or $y$), or 
reaction with different sets of atoms on right and left sides of 
the equation.

Thus, the list of 9181 valid reactions was formed. 
Then automatic balancing and calculations of the masses were performed 
for each reaction on this list. The results (benched on Ubuntu 24.04, 
AMD Ryzen 7 5700x, and 64 GB DDR4 RAM with output writing to a .txt file) are 0.057 ms
per formula and 0.49 ms per reaction. There are 69 million reactions in the Reaxys 
database [@ReaxysDatabase], therefore all the reactions in the largest database can be balanced within
9.5 hours on a single consumer-grade PC!

# References