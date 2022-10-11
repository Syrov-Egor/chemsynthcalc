ChemicalFormula class
=====================
A base class to represent a single chemical formula.

:class:`chemsynthcalc.chemical_formula.ChemicalFormula`

Formula object initialization
-----------------------------
To create a formula object from a string::

    from chemsynthcalc import ChemicalFormula

    formula_string = "(NH4)2SO4"
    formula = ChemicalFormula(formula_string)

.. important::
    The allowed symbols for ChemicalFormula input string are:
    
    a-z A-Z 0-9 . () {} [] * · • *whitespace*
    
    Whitespaces will be ignored. If there are any other symbols 
    in the string, they won't be ignored, instead
    :class:`chemsynthcalc.chem_errors.InvalidCharacter` exception will be raised.

There is an optional **rounding_order (int)** parameter for rounding precision::
    
    from chemsynthcalc import ChemicalFormula

    formula_string = "(NH4)2SO4"
    formula = ChemicalFormula(
        formula = formula_string,
        rounding_order = 8
        )

Formula object properties
-------------------------
After the object initialization, we can access ChemicalFormula properties:

* parsed_formula
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.parsed_formula`

* molar_mass 
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.molar_mass`

* mass_percent
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.mass_percent`

* atomic_percent
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.atomic_percent`

* oxide_percent
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.oxide_percent`

* output_results
    :py:attr:`chemsynthcalc.chemical_formula.ChemicalFormula.output_results`

Output
-------------------------
One can output ChemicalFormula results using one of the 4 methods:

* print_results: print to console
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.print_results()`

* export_to_txt: save as plain txt file
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.export_to_txt()`

* as_json: serialization of output into JSON object
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.as_json()`

* export_to_json: save as JSON file
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.export_to_json()`