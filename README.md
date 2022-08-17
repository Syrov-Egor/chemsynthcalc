# Chemical synthesis calculator
This is a Python package for the calculation of the necessary masses of substances for chemical synthesis directly from the reaction string. It includes solutions for all intermidiate steps, including chemical formula parsing, molar mass calculation and reaction balancing with different matrix methods.

## Prerequisites
* [Python](https://www.python.org/downloads/) 3.7+
* [NumPy](https://numpy.org/)
* [SciPy](https://scipy.org/)

## Installation
Install from [pypi](https://pypi.org/):

`pip install chemsynthcalc`

## Usage
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

