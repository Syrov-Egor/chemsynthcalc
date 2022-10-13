ChemicalReaction class
======================
An interface class to parse, decompose and calculate end results of
chemical reaction from the reaction string.

:class:`chemsynthcalc.chemical_reaction.ChemicalReaction`


Reaction object initialization
------------------------------
To create a reaction object from a string::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string)

.. important::
    The symbols allowed for the ChemicalReaction input string are:
    
    a-z A-Z 0-9 . () [] {} * · • = < - > → ⇄ + *whitespace*
    
    Whitespaces will be ignored. If there are any other symbols 
    in the string, they will not be ignored, instead
    :class:`chemsynthcalc.chem_errors.InvalidCharacter` exception will be raised.

.. important::
    Other important conditions for calculations are:

    There must be a reactant-product separator (listed in :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.possible_reaction_separators`)
    
    The sets of atoms in the reactant and product (left-right) parts of equation should be equal 
    (in other words, all atoms in the left must appear on the right). If not, :class:`chemsynthcalc.chem_errors.ReactantProductDifference` \
    exception will be raised.

There are other important arguments for reaction object creation. The most important
concept is calculation *modes*.

Modes
-----

force mode
++++++++++

.. warning::
    This mode does not imply any automatic checking of calculation results.
    Use it at your own risk!

Force mode can be used when user enters coefficients in the reaction string
and want masses to be calculated whether the reaction is balanced or not.
For example, consider the following reaction for iodine synthesis::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"

If we know beforehand that coefficient list will be [1,5,3,3,3,3], we can input::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="force")

We can make sure that the reaction is balanced::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="force")
    print(reaction.is_balanced)
    # True

And calculate the respective masses from this reaction::
        
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="force")
    print(reaction.is_balanced)
    print(reaction.masses)

    # True
    # [0.28105463, 1.09008406, 0.3864145, 1.0, 0.6865721, 0.07098109]

But what if we decided to add a 3-fold excess of KIO3? In the force mode, we **can**
do that and still get output masses, even though the reaction is not balanced::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "3KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="force")
    print(reaction.is_balanced)
    print(reaction.masses)

    # False
    # [0.84316391, 1.09008406, 0.3864145, 1.0, 0.6865721, 0.07098109]

As we can see, the mass of KIO3 has increased, while all other masses
have not been changed.

check mode
++++++++++
Check mode is almost the same as force mode, but with mandatory reaction balance checks.
Therefore, an unbalanced reaction cannot be calculated::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "3KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="check")
    print(reaction.is_balanced)
    print(reaction.masses)
    
    # chemsynthcalc.chem_errors.ReactionNotBalanced: This reaction is not balanced!

Whereas if we input the right list of coefficients::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
    reaction = ChemicalReaction(reaction_string, mode="check")
    print(reaction.is_balanced)
    print(reaction.masses)
    
    # True
    # [0.28105463, 1.09008406, 0.3864145, 1.0, 0.6865721, 0.07098109]

balance mode
++++++++++++
The last mode is  designed for the automatic reaction balancing by
means of :class:`chemsynthcalc.reaction_balance.Balancer` class.
For most reactions, the auto balancing option is enough::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    print(reaction.coefficients)
    print(reaction.is_balanced)
    print(reaction.algorithm) # we can also check which algorithm solved the reaction
    print(reaction.masses)
    
    # [1, 5, 3, 3, 3, 3]
    # True
    # inverse
    # [0.28105463, 1.09008406, 0.3864145, 1.0, 0.6865721, 0.07098109]

In some cases, however, the auto-balancing is not enough, or one would
want to calculate coefficients strictly with a specific algorithm. To 
address these issues, the following is implemented in ChemicalReaction class
logic:

Coefficients property setter
----------------------------
Unlike other properties of ChemicalFormula and ChemicalReaction, which are
programmed as get-only cached properties without setters,
the :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients` property
can be set directly with appropriate list of coefficients::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
    reaction = ChemicalReaction(reaction_string, mode="check")
    reaction.coefficients = [1, 5, 3, 3, 3, 3]
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # [1, 5, 3, 3, 3, 3]
    # True

The :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients_check()` method checks
setted coefficients *after* setting in such properties as 
:py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.normalized_coefficients` and
:py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.masses`, therefore does not allow
methods to calculate values with bad coefficients, for example::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
    reaction = ChemicalReaction(reaction_string, mode="check")
    reaction.coefficients = [1]
    print(reaction.coefficients)
    print(reaction.is_balanced)
    print(reaction.masses)

    # [1]
    # False
    # chemsynthcalc.chem_errors.BadCoeffiecients: number of coefficients should be equal to 6

This is one more way (along with direct input of coefficients in the reaction string for force
and check modes) to input a custom set of coefficients for a particular reaction.

.. note::
    Coefficients setting works in all three modes (force, check and balance).

Coefficients property calculation
---------------------------------
The :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.balance_reaction()` method
is designed to give high-level interface to every coefficient calculation algorithm
implemented in chemsynthcalc. These can be chosen by the first argument, "algorithm",
and they are:

inv or matrix inverse Thorne algorithm
++++++++++++++++++++++++++++++++++++++
See :meth:`chemsynthcalc.reaction_balance.Balancer.inv_algorithm()`

gpinv or general pseudoinverse Risteski algorithm
+++++++++++++++++++++++++++++++++++++++++++++++++
See :meth:`chemsynthcalc.reaction_balance.Balancer.gpinv_algorithm()`

ppinv or partial pseudoinverse Risteski algorithm
+++++++++++++++++++++++++++++++++++++++++++++++++
See :meth:`chemsynthcalc.reaction_balance.Balancer.ppinv_algorithm()`

comb or combinatorial algorithm
+++++++++++++++++++++++++++++++++++++++++++++++++
See :meth:`chemsynthcalc.reaction_balance.Balancer.comb_algorithm()`

One can also specify whether or not calculated coefficients should be integers
(by "intify" bool flag), and the maximum number of combinations ("max_comb") for 
the combinatorial algorithm.

.. note::
    Although :meth:`chemsynthcalc.reaction_balance.Balancer.intify_coefficients()`
    *tries* to intify coefficients, it won't always succeed. In this case, coefficients
    will stay *float* without any exception raise.

Some examples of the same reaction balanced by these four different methods::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.coefficients = reaction.balance_reaction(algorithm='inv')
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # Can't equalize this reaction by inverse algorithm
    # None
    # False

The inverse algorithm can't handle this! Let's try others, like general pseudoinverse::
        
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.coefficients = reaction.balance_reaction(algorithm='gpinv')
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # [1954, 1854, 518, 1093, 1096, 55, 901, 898]
    # True

partial pseudoinverse::
            
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.coefficients = reaction.balance_reaction(algorithm='ppinv')
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # [39, 39, 13, 13, 11, 5, 16, 18]
    # True

and combinatorial::
       
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.coefficients = reaction.balance_reaction(algorithm='comb')
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # [4, 5, 1, 1, 1, 1, 1, 3]
    # True

As we can see, we have got *four* different results (including 3 right ones) using *four* different algorithms.
This is why the :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.balance_reaction()`
method was implemented in the first place. We can, of course, get gpinv or ppinv data without intification::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.coefficients = reaction.balance_reaction(algorithm='gpinv', intify=False)
    print(reaction.coefficients)
    print(reaction.is_balanced)

    # [1.4169688179840456, 1.3444525018129083, 0.37563451776649776, 0.7926033357505439, 0.7947788252356772, 0.03988397389412329, 0.6533720087019586, 0.6511965192168249]
    # True

.. note::
    Coefficients calculation by :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.balance_reaction()` 
    works only in the balance mode.

Now, when we covered *modes* and coefficients property calculations, we can get back to other
arguments for reaction object creation.

Target
------
The target is the target chemical synthesis substance (i.e., whose mass is known in advance).
The target choice is implemented as an integer pointer to formula index in products (right side of equation).
There are, of course, other ways to do this (like explicitly input target as formula string),
but integer index target was chose as method less prone to errors. Most of the time the target
of the synthesis is the first product anyway (which is equal to 0 by default).
We can set a target with initialization (the target is 1 which is FeO)::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance", target=1)
    print(reaction.masses)

    # [3.97359366, 0.28358172, 1.52731369, 1.0, 0.77944268, 0.12575572, 0.32138621, 0.5032771]

If we change the target, masses will obviously change too::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance", target=2)
    print(reaction.masses)

    # [5.09799345, 0.36382626, 1.95949455, 1.28296797, 1.0, 0.16134056, 0.41232821, 0.64568841]

Target mass
-----------
The mass of target compound, the only mass that we know in advance before synthesis (in grams).
1 gram by default. We can change the target mass during the ChemicalReaction object
initialization::

    from chemsynthcalc import ChemicalReaction

    reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
    reaction = ChemicalReaction(reaction_string, mode="balance", target=1, target_mass=2)

    print(reaction.masses)
    
    # [7.94718732, 0.56716344, 3.05462737, 2.0, 1.55888536, 0.25151145, 0.64277241, 1.0065542]

As we can see, mass of the target (FeO) is now 2.0, and all other masses also have been changed.

Other arguments
---------------
There are two more arguments of ChemicalReaction object.

The *try_comb* bool flag determines whether there should be an attempt to calculate the coefficients
using the combinatorial method in automatic balance mode if all other methods have failed.

*rounding_order* (int) parameter for rounding precision.

ChemicalReaction properties
---------------------------
After the object initialization, we can access the ChemicalReaction properties:

* reaction
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.reaction`

* separator
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.separator`

* reactants
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.reactants`

* products
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.products`

* compounds
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.compounds`

* initial_coefficients
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.initial_coefficients`

* formulas
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.formulas`

* matrix
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.matrix`

* reactant_matrix
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.reactant_matrix`

* product_matrix
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.product_matrix`

* molar_masses
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.molar_masses`

* coefficients
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients`

* normalized_coefficients
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.normalized_coefficients`

* is_balanced
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.is_balanced`

* final_reaction
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.final_reaction`

* final_reaction_normalized
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.final_reaction_normalized`

* masses
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.masses`

* output_results
    :py:attr:`chemsynthcalc.chemical_reaction.ChemicalReaction.output_results`

Output
------
A typical ChemicalReaction results output will look like this::
    
    from chemsynthcalc import ChemicalReaction

    reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
    reaction = ChemicalReaction(reaction_string, mode="balance")
    reaction.print_results()

    initial reaction: KIO3+KI+H2SO4=I2+K2SO4+H2O
    reaction matrix:
    [[1. 1. 0. 0. 2. 0.]
     [1. 1. 0. 2. 0. 0.]
     [3. 0. 4. 0. 4. 1.]
     [0. 0. 2. 0. 0. 2.]
     [0. 0. 1. 0. 1. 0.]]
    mode: balance
    coefficients: [1, 5, 3, 3, 3, 3]
    normalized coefficients: [0.33333333, 1.66666667, 1, 1, 1, 1]
    balanced by algorithm: inverse
    is balanced: True
    final reaction: KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O
    final reaction normalized: 0.33333333KIO3+1.66666667KI+H2SO4=I2+K2SO4+H2O
    target: I2
    KIO3: M = 213.9950 g/mol, m = 0.2811 g
    KI: M = 165.9980 g/mol, m = 1.0901 g
    H2SO4: M = 98.0720 g/mol, m = 0.3864 g
    I2: M = 253.8000 g/mol, m = 1.0000 g
    K2SO4: M = 174.2520 g/mol, m = 0.6866 g
    H2O: M = 18.0150 g/mol, m = 0.0710 g

One can output ChemicalReaction results using one of the 4 methods:

* print_results: print to console
    :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.print_results()`

* export_to_txt: save as plain txt file
    :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.export_to_txt()`

* as_json: serialization of output into an JSON object
    :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.as_json()`

* export_to_json: save as an JSON file
    :meth:`chemsynthcalc.chemical_reaction.ChemicalReaction.export_to_json()`