# Chemical synthesis calculator
A Python package for the calculation of the necessary masses of substances for chemical synthesis directly from the reaction string. It includes solutions for all intermidiate steps, including chemical formula parsing, molar mass calculation and reaction balancing with different matrix methods.

## Prerequisites
* [Python](https://www.python.org/downloads/) 3.8+
* [NumPy](https://numpy.org/)

## Installation
Install from [pypi](https://pypi.org/):

`pip install chemsynthcalc`

## Example use
Let's say that we need to prepare 3 grams of [YBCO](https://en.wikipedia.org/wiki/Yttrium_barium_copper_oxide) by solid-state synthesis from respective carbonates. The reaction string will look something like this (to simplify, let's leave it without oxygen nonstoichiometry):

```Python
reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"
```

Now, we can create a chemical reaction object of `ChemicalReaction` class, which will be used in the calculation. We need to specify arguments for our particular case:
```Python
from chemsynthcalc import ChemicalReaction

reaction = ChemicalReaction(
    reaction = reaction_string, # our reaction string
    target = 0, # number of target compound in product list
    target_mass = 3 # desired mass of target compound,
    mode = "balance" # mode of coefficients calculations,
)
```

Now, to perform automatic calculation, all we need to do is to put:
```Python
reaction.print_results(print_rounding_order=4) # assuming we use analytical balances with 4 digits presicion
```

And we get our ouput in the console:
```
initial reaction: BaCO3+Y2(CO3)3+CuCO3+O2→YBa2Cu3O7+CO2
reaction matrix:
[[1. 0. 0. 0. 2. 0.]
 [0. 2. 0. 0. 1. 0.]
 [0. 0. 1. 0. 3. 0.]
 [1. 3. 1. 0. 0. 1.]
 [3. 9. 3. 2. 7. 2.]]
mode: balance
coefficients: [8, 2, 12, 1, 4, 26]
normalized coefficients: [2, 0.5, 3, 0.25, 1, 6.5]
balanced by Thorne algorithm
target: YBa2Cu3O7
final reaction: 8BaCO3+2Y2(CO3)3+12CuCO3+O2→4YBa2Cu3O7+26CO2
final reaction normalized: 2BaCO3+0.5Y2(CO3)3+3CuCO3+0.25O2→YBa2Cu3O7+6.5CO2
BaCO3: M = 197.3380 g/mol, m = 1.7773 g
Y2(CO3)3: M = 357.8360 g/mol, m = 0.8057 g
CuCO3: M = 123.5540 g/mol, m = 1.6692 g
O2: M = 31.9980 g/mol, m = 0.0360 g
YBa2Cu3O7: M = 666.1970 g/mol, m = 3.0000 g
CO2: M = 44.0090 g/mol, m = 1.2882 g
```
Now, we got all masses ready for our planned synthesis!

See the docmentation for detailed functionality and all features of the package.

## Functionality
### ChemicalFormula class
To create a single chemical formula object: 
```Python
from chemsynthcalc import ChemicalFormula

formula = ChemicalFormula("K2SO4")
```
Now, we can get properties of this object, which are:

* `parsed_formula`  
  A dictionary machine-readable representation of chemical formula needed for futher calculations. The keys are atoms in formula and the values are their coefficients. **The resulting dictionary is unordered**.
  ```Python
  >>> formula.parsed_formula
  {'K': 2.0, 'O': 4.0, 'S': 1.0}
  ```
* `molar_mass`  
  A float number of [molar mass](https://en.wikipedia.org/wiki/Molar_mass) of the compound in g/mol.
  ```Python
  >>> formula.molar_mass
  174.2592
  ```
* `mass_percent`  
  A dictionary of [relative mass fraction](https://en.wikipedia.org/wiki/Mass_fraction_(chemistry)) of each atom with 100% sum.
  ```Python
  >>> formula.mass_percent
  {'S': 18.40075, 'K': 44.87373, 'O': 36.72552}
  ```
* `atomic_percent`  
  A dictionary of [relative mole fraction](https://en.wikipedia.org/wiki/Mole_fraction) of each atom with 100% sum.
  ```Python
  >>> formula.atomic_percent
  {'S': 14.28571, 'K': 28.57143, 'O': 57.14286}
  ```
* `oxide_percent`  
  A dictionary of [oxide fraction](https://d32ogoqmya1dw8.cloudfront.net/files/introgeo/studio/examples/minex02.pdf) of metals with 100% sum. Oxide types are specified in `periodic_table.py`.
  ```Python
  >>> formula.oxide_percent
  {'K2O': 54.05511, 'SO3': 45.94489}
  ```
### ChemicalReaction class
To create an object that represents a chemical reaction:
```Python
from chemsynthcalc import ChemicalReaction

reaction_string = "KMnO4+MnSO4+H2O=MnO2+K2SO4+H2SO4"
reaction = ChemicalReaction(reaction_string)
```
## License
The code is provided under the MIT license.

