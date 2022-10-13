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
    The symbols allowed for the ChemicalFormula input string are:
    
    a-z A-Z 0-9 . () {} [] * · • *whitespace*
    
    Whitespaces will be ignored. If there are any other symbols 
    in the string, they will not be ignored, instead
    :class:`chemsynthcalc.chem_errors.InvalidCharacter` exception will be raised.
    
    All brackets in formula should be paired bracket-type-wise (), {}, [].
    If not, :class:`chemsynthcalc.chem_errors.BracketsNotPaired` exception will be raised.

.. important::
    Only one adduct (like X*H2O) per formula is allowed. If parser detects more
    than one adduct symbols (\*·•) it will raise :class:`chemsynthcalc.chem_errors.MoreThanOneAdduct`
    exception.

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
A typical ChemicalFormula results output will look like this::
    
    from chemsynthcalc import ChemicalFormula

    ChemicalFormula("(NH4)2SO4").print_results()

    formula: (NH4)2SO4
    parsed formula: {'N': 2.0, 'H': 8.0, 'S': 1.0, 'O': 4.0}
    molar mass: 132.134
    mass percent: {'N': 21.2012, 'H': 6.1029, 'S': 24.2632, 'O': 48.4327}
    atomic percent: {'N': 13.3333, 'H': 53.3333, 'S': 6.6667, 'O': 26.6667}
    oxide percent: {'NO2': 37.6894, 'H2O': 29.5174, 'SO3': 32.7932}

One can output ChemicalFormula results using one of the 4 methods:

* print_results: print to console
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.print_results()`

* export_to_txt: save as plain txt file
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.export_to_txt()`

* as_json: serialization of output into an JSON object
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.as_json()`

* export_to_json: save as an JSON file
    :meth:`chemsynthcalc.chemical_formula.ChemicalFormula.export_to_json()`