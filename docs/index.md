#About

**ChemSynthCalc** stands for **Chem**ical **Synth**esis **Calc**ulator - 
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

##Main project page

[https://github.com/Syrov-Egor/chemsynthcalc](https://github.com/Syrov-Egor/chemsynthcalc)

##Background

ChemSynthCalc was created to address the issue of inorganic synthesis calculations.
There are three main aspects to this problem. 

First of all, a single student, scientist, postdoc, etc., sitting in their lab 
trying to synthesize a new material. The sample size of potential reactions is quite small in this case, 
and our fellow scientist can use the old-fashioned pen and paper method. 
But why should one not try to automate this boring task?

Secondly, while we are not quite there yet, one can imagine a robotic inorganic synthesis
station, like [Dr. Cronin](https://pubs.rsc.org/en/content/articlelanding/2012/LC/c2lc40761b)
devices. In this case, we need to calculate a large number of reactions fast and precisely, 
and we really should not hard-code all the masses beforehand.

Finally, there are datasets of [text-mined inorganic reactions](https://www.nature.com/articles/s41597-019-0224-1).
They will surely expand and grow in size. With a sample size of tens of thousands, we need our reaction balancing software 
to be extremely fast, robust, and flexible enough to balance as many reactions as it can. 
In this case, we surely need a free open-source solution that can be embedded in the data processing system.

ChemSynthCalc addresses all three of those cases. It is simple enough to use for a single scientist who 
is familiar with using Python packages and fast and robust enough to precisely 
calculate hundreds and thousands of reactions.

##Competitive analysis

There are already a large number of reaction balancing software, why do you need
chemsynthcalc and why is it better than competitors?

* chemsynthcalc is completely free and open-source (under MIT license)
* chemsynthcalc provides a rich and simple API for its functions
* chemsynthcalc can balance a huge variety of reactions
* chemsynthcalc can deal with formulas with float atom count (like RbLa~0.99~Eu~0.01~Nb~2~O~7~), both for molar mass and reaction balancing
* chemsynthcalc supports an infinite amount of nested parentheses in formulas
* chemsynthcalc supports adduct notaion (like CuSO~4~*5H~2~O)
* chemsynthcalc is fast thanks to NumPy matrix operations
* finally, chemsynthcalc is one of the few programs that can directly output precursor masses from a reaction string (in three lines of code!)

|Software                                                    |Mass calculation|Float coefficients and atom amounts|API|
|------------------------------------------------------------|:----:|:----:|:----:|
|[chemsynthcalc](https://github.com/Syrov-Egor/chemsynthcalc)|<span style="color:green">✓</span>|<span style="color:green">✓</span>|<span style="color:green">✓</span>|
|[WebQC](https://www.webqc.org/balance.php)|<span style="color:green">✓</span>|<span style="color:red">x</span>|<span style="color:red">x</span>|
|[Equation Balancer](https://equationbalancer.com/)|<span style="color:red">x</span>|<span style="color:red">x</span>|<span style="color:red">x</span>|
|[ChemicalAid](https://en.intl.chemicalaid.com/tools/equationbalancer.php)|<span style="color:red">x</span>|<span style="color:green">✓</span>|<span style="color:red">x</span>|
|[CHEMIX School](https://www.chemix-chemistry-software.com/chemistry-software.html)|<span style="color:red">x</span>|<span style="color:green">✓</span>|<span style="color:red">x</span>|
|[ChemPy](https://github.com/bjodah/chempy)|<span style="color:red">x</span>|<span style="color:red">x</span>|<span style="color:green">✓</span>|

There are two main classes to acesses the functionality: ChemicalFormula and ChemicalReaction. You can check them in the menu.