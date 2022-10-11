ChemicalReaction class
======================
An interface class to parse, decompose and calculate end results of
chemical reaction from the reaction string.

:class:`chemsynthcalc.chemical_reaction.ChemicalReaction`


Object initialization
---------------------
To create a reaction object from a string::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string)

.. important::
    The allowed symbols for ChemicalReaction input string are:
    
    a-z A-Z 0-9 . () [] {} * · • = < - > → ⇄ + *whitespace*
    
    Whitespaces will be ignored. If there are any other symbols 
    in the string, they won't be ignored, instead
    :class:`chemsynthcalc.chem_errors.InvalidCharacter` exception will be raised.

There are others important arguments for reaction object creation. The most important
concept is calculation *modes*.

Modes
-----

force mode
++++++++++

check mode
++++++++++

balance mode
++++++++++++