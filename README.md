# Chemical synthesis calculator
Python package for calculating the masses of substances required for chemical synthesis directly from the reaction string. It includes solutions for all intermediate steps, including chemical formula parsing, molar mass calculation and reaction balancing with different matrix methods.

## Prerequisites
* [Python](https://www.python.org/downloads/) 3.10+
* [SciPy](https://scipy.org/)

## Installation
Install from [pypi](https://pypi.org/project/chemsynthcalc/):

``` Bash
pip install chemsynthcalc
```

## Documentation
[Read the documentation](https://syrov-egor.github.io/chemsynthcalc/)

## Example use
Let's say that we need to prepare 3 grams of [YBCO](https://en.wikipedia.org/wiki/Yttrium_barium_copper_oxide) by solid-state synthesis from respective carbonates. The reaction string will look something like this (to simplify, let's leave it without oxygen nonstoichiometry):

``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"
```

Now, we can create a chemical reaction object of the `ChemicalReaction` class, which will be used in the calculation. We need to specify the arguments for our particular case:
``` Python
>>> reaction = ChemicalReaction(
    reaction = reaction_string, # our reaction string
    target = 0, # index of target compound in the product list
    target_mass = 3, # desired mass of target compound,
    mode = "balance" # mode of coefficients calculations,
)
```

Now, to perform the automatic calculation, all we need to do is to put:
``` Python
>>> reaction.print_results(print_rounding_order=4) 
# assuming that we use analytical balances with 4 digit presicion
```

And we get our output in the terminal:
``` Python
initial reaction: BaCO3+Y2(CO3)3+CuCO3+O2→YBa2Cu3O7+CO2
reaction matrix:
 [[1. 0. 0. 0. 2. 0.]
 [1. 3. 1. 0. 0. 1.]
 [3. 9. 3. 2. 7. 2.]
 [0. 2. 0. 0. 1. 0.]
 [0. 0. 1. 0. 3. 0.]]
mode: balance
formulas: ['BaCO3', 'Y2(CO3)3', 'CuCO3', 'O2', 'YBa2Cu3O7', 'CO2']
coefficients: [8, 2, 12, 1, 4, 26]
normalized coefficients: [2, 0.5, 3, 0.25, 1, 6.5]
algorithm: inverse
is balanced: True
final reaction: 8BaCO3+2Y2(CO3)3+12CuCO3+O2→4YBa2Cu3O7+26CO2
final reaction normalized: 2BaCO3+0.5Y2(CO3)3+3CuCO3+0.25O2→YBa2Cu3O7+6.5CO2
molar masses: [197.335, 357.835676, 123.554, 31.998, 666.190838, 44.009]
target: YBa2Cu3O7
masses: [1.7773, 0.8057, 1.6692, 0.036, 3.0, 1.2882]
BaCO3: M = 197.3350 g/mol, m = 1.7773 g
Y2(CO3)3: M = 357.8357 g/mol, m = 0.8057 g
CuCO3: M = 123.5540 g/mol, m = 1.6692 g
O2: M = 31.9980 g/mol, m = 0.0360 g
YBa2Cu3O7: M = 666.1908 g/mol, m = 3.0000 g
CO2: M = 44.0090 g/mol, m = 1.2882 g
```
Now, we've got all the masses ready for our planned synthesis!

## Features
* Formula parsing
  ```Python
  >>> from chemsynthcalc import ChemicalFormula

  >>> ChemicalFormula("C2H5OH").parsed_formula
  {'C': 2.0, 'H': 6.0, 'O': 1.0}
  ```
* Calculation of the molar mass 
  ```Python
  >>> ChemicalFormula("C2H5OH").molar_mass
  46.069
  ```
* [Mass](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)), [atomic](https://en.wikipedia.org/wiki/Mole_fraction), and [oxide](https://d32ogoqmya1dw8.cloudfront.net/files/introgeo/studio/examples/minex02.pdf) percent calculations with `mass_percent`, `atomic_percent` and `oxide_percent` (including user-defined oxides) properties of `ChemicalFormula`.
* Auto-balancing chemical equations by 4 different matrix methods in `balance` mode:
  ```Python
  >>> from chemsynthcalc import ChemicalReaction

  >>> reaction_string = "K4Fe(CN)6 + KMnO4 + H2SO4 = KHSO4 + Fe2(SO4)3 + MnSO4 + HNO3 + CO2 + H2O"

  >>> ChemicalReaction(reaction_string, mode="balance").final_reaction
  10K4Fe(CN)6+122KMnO4+299H2SO4=162KHSO4+5Fe2(SO4)3+122MnSO4+60HNO3+60CO2+188H2O
  ```
* Calculation of masses for user-defined coefficients in `force` (calculates regardless of balance) and `check` (checks if reaction is balanced by user-defined coefficients) modes.
  ```Python
  >>> ChemicalReaction("BaCO3+TiO2=BaTiO3", mode="force").masses #we can drop CO2 product and still get masses in this mode. 
  [0.84623763, 0.34248749, 1.0]
  ```
  ```Python
  >>> ChemicalReaction("2H2+O2=2H2O", mode="check").coefficients #we can be sure that reaction is balanced with our coefficients in this mode
  [2, 1, 2]
  ```
  Setting the coefficients directly into `ChemicalReaction` instance:
  ```Python
  >>> reaction = ChemicalReaction("H2+O2=H2O", mode="check")
  >>> reaction.coefficients = [2,1,2]
  >>> reaction.coefficients
  [2, 1, 2]
  >>> reaction.is_balanced
  True
  ```
* Calculation of coefficients with `ChemicalReaction.balance` object individually by each of 4 different algorithms (inverse, general pseudoinverse, partial pseudoinverse and combinatorial algorithms).
* Export of results of both `ChemicalFormula` and `ChemicalReaction` into .txt file (with `.to_txt()`), into JSON object (with `.to_json()`) or JSON file (with `.to_json_file()`).

## License
The code is provided under the MIT license.

## Contact
If you have any questions, please contact **Egor Syrov** at syrov_ev@mail.ru or
create an issue at github https://github.com/Syrov-Egor/chemsynthcalc/issues.