A class to represent a single chemical formula.

API docs: [ChemicalFormula][chemsynthcalc.chemical_formula.ChemicalFormula]

## Formula object initialization

To create a formula object from a string:

``` Python
from chemsynthcalc import ChemicalFormula

formula_string = "(NH4)2SO4"
formula = ChemicalFormula(formula_string)
```

!!! warning "Important"

    The symbols allowed for the ChemicalFormula input string are:
    a-z A-Z 0-9 . () {} [] * · • whitespace

Whitespaces will be ignored. If there are any other symbols in the string, they will not be ignored, instead [InvalidCharacter][chemsynthcalc.chem_errors.InvalidCharacter] error will be raised.

All brackets in formula should be paired bracket-type-wise (), {}, []. If not, [BracketsNotPaired][chemsynthcalc.chem_errors.BracketsNotPaired] error will be raised.

!!! warning "Important"

    Only one adduct (like X\*H2O) per formula is allowed. If parser detects more than one adduct symbols (*·•) it will raise [MoreThanOneAdduct][chemsynthcalc.chem_errors.MoreThanOneAdduct] error.

There is an optional rounding_order (int) parameter for rounding precision:

``` Python
from chemsynthcalc import ChemicalFormula

formula_string = "(NH4)2SO4"
formula = ChemicalFormula(
    formula = formula_string,
    rounding_order = 8
    )
```

## Formula object properties

After the object initialization, we can access ChemicalFormula properties:

* [formula][chemsynthcalc.chemical_formula.ChemicalFormula]
* [parsed_formula][chemsynthcalc.chemical_formula.ChemicalFormula.parsed_formula]
* [molar_mass][chemsynthcalc.chemical_formula.ChemicalFormula.molar_mass]
* [mass_percent][chemsynthcalc.chemical_formula.ChemicalFormula.mass_percent]
* [atomic_percent][chemsynthcalc.chemical_formula.ChemicalFormula.atomic_percent]
* [oxide_percent][chemsynthcalc.chemical_formula.ChemicalFormula.oxide_percent]
* [output_results][chemsynthcalc.chemical_formula.ChemicalFormula.output_results]

## Custom oxides

The [oxide_percent][chemsynthcalc.chemical_formula.ChemicalFormula.oxide_percent] property calculates relative percentages of oxides of 
non-oxygen elements in the formula (tipycally metals). The default oxide formulas are listed in [PeriodicTable][chemsynthcalc.periodic_table.PeriodicTable]. One can feed their custom oxide formulas while instantiating an object:

``` Python
>>> from chemsynthcalc import ChemicalFormula

>>> ChemicalFormula("(NH4)2SO4").oxide_percent
{'NO2': 37.68939937, 'H2O': 29.51742331, 'SO3': 32.79317732}

>>> ChemicalFormula("(NH4)2SO4", "S2O3").oxide_percent
{'NO2': 41.79831326, 'H2O': 32.73542499, 'S2O3': 25.46626175}
```

## Output

A typical ChemicalFormula results output will look like this:

``` Python
>>> from chemsynthcalc import ChemicalFormula

>>> ChemicalFormula("(NH4)2SO4").print_results()

formula: (NH4)2SO4
parsed formula: {'N': 2.0, 'H': 8.0, 'S': 1.0, 'O': 4.0}
molar mass: 132.134
mass percent: {'N': 21.2012, 'H': 6.1029, 'S': 24.2632, 'O': 48.4327}
atomic percent: {'N': 13.3333, 'H': 53.3333, 'S': 6.6667, 'O': 26.6667}
oxide percent: {'NO2': 37.6894, 'H2O': 29.5174, 'SO3': 32.7932}
```

One can output ChemicalFormula results using one of the 4 methods:

* [print_results][chemsynthcalc.chemical_formula.ChemicalFormula.print_results]: print to stdout
* [to_txt][chemsynthcalc.chemical_formula.ChemicalFormula.to_txt]: save as plain txt file
* [to_json][chemsynthcalc.chemical_formula.ChemicalFormula.to_json]: serialization of output into an JSON object
* [to_json_file][chemsynthcalc.chemical_formula.ChemicalFormula.to_json_file]: save as JSON file