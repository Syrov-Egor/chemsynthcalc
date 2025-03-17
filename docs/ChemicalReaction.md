A class to calculate end results of chemical reaction from the reaction string.

API docs: [ChemicalReaction][chemsynthcalc.chemical_reaction.ChemicalReaction]

## Reaction object initialization

To create a reaction object from a string:

``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
reaction = ChemicalReaction(reaction_string)
```

!!! warning "Important"

    The symbols allowed for the ChemicalReaction input string are:
    a-z A-Z 0-9 . () [] {} * · • = < - > → ⇄ + whitespace

Whitespaces will be ignored. If there are any other symbols in the string, they will not be ignored, instead [InvalidCharacter][chemsynthcalc.chem_errors.InvalidCharacter] error will be raised.

!!! warning "Important"
    Other important conditions for calculations are:

    There must be one reactant-product separator (listed in [possible_reaction_separators][chemsynthcalc.reaction.Reaction])

    The sets of atoms in the reactant and product (left-right) parts of equation should be equal (in other words, all atoms in the left must appear on the right). If not, [ReactantProductDifference][chemsynthcalc.chem_errors.ReactantProductDifference] error will be raised.

There are other important arguments for reaction object creation. The most important concept is calculation *mode*.

## Modes

### Force mode

!!! danger

    This mode does not imply any automatic checking of calculation results. Use it at your own risk!

Force mode can be used when user enters coefficients in the reaction string and want masses to be calculated whether the reaction is balanced or not. For example, consider the following reaction for iodine synthesis:

``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
```

If we know beforehand that coefficient list will be [1,5,3,3,3,3], we can input:

``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
reaction = ChemicalReaction(reaction_string, mode="force")
```
We can make sure that the reaction is balanced:

``` Python
>>> reaction.is_balanced
True
```
And calculate the respective masses from this reaction:

``` Python
>>> reaction.masses
[0.2810506, 1.09007501, 0.38640089, 1.0, 0.68654792, 0.07097859]
```

But what if we decided to add a 3-fold excess of KIO3? In the force mode, we **can** do that and still get output masses, even though the reaction is not balanced:

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "3KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="force")

>>> reaction.is_balanced
False

>>> reaction.masses
[0.84315182, 1.09007501, 0.38640089, 1.0, 0.68654792, 0.07097859]
```

As we can see, the mass of KIO3 has increased, while all other masses have not been changed.

### Check mode

Check mode is almost the same as force mode, but with mandatory reaction balance checks. Therefore, an unbalanced reaction cannot be calculated:

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "3KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="check")
>>> reaction.is_balanced

chemsynthcalc.chem_errors.ReactionNotBalanced: This reaction is not balanced!
```

Whereas if we input the right list of coefficients:
``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "3KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="check")

>>> reaction.is_balanced
True

>>> reaction.masses
[0.2810506, 1.09007501, 0.38640089, 1.0, 0.68654792, 0.07097859]
```

### Balance mode
The last mode is designed for the automatic reaction balancing by means of [Balancer][chemsynthcalc.balancer.Balancer] class. It is the **default mode**. For most reactions, the auto balancing option is enough:


``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="balance")

>>> reaction.coefficients
[1, 5, 3, 3, 3, 3]

>>> reaction.is_balanced
True

>>> reaction.algorithm # we can also check which algorithm solved the reaction
inverse

>>> reaction.masses
[0.2810506, 1.09007501, 0.38640089, 1.0, 0.68654792, 0.07097859]
```

In some cases, however, the auto-balancing is not enough, or one would want to calculate coefficients strictly with a specific algorithm. To address these issues, the following is implemented in ChemicalReaction class logic:

## Coefficients property calculation

The [Balancer][chemsynthcalc.balancer.Balancer] class is designed to give high-level interface to every coefficient calculation algorithm implemented in chemsynthcalc. These can be chosen by the specific method name, and they are:

### inv or matrix inverse Thorne algorithm
See [inv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._inv_algorithm] for details.

### gpinv or general pseudoinverse Risteski algorithm
See [gpinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._gpinv_algorithm] for details.

### ppinv or partial pseudoinverse Risteski algorithm
See [ppinv_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._ppinv_algorithm] for details.

### comb or combinatorial algorithm
See [comb_algorithm][chemsynthcalc.balancing_algos.BalancingAlgorithms._comb_algorithm] for details.


!!! note

    Although [_intify_coefficients][chemsynthcalc.balancer.Balancer._intify_coefficients] tries to intify coefficients, it won't always succeed. In this case, coefficients will stay float without any exception raise.

Some examples of the same reaction balanced by these four different methods:

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
>>> reaction = ChemicalReaction(reaction_string, mode="balance")

>>> reaction.coefficients = reaction.balancer.inv()
"chemsynthcalc.chem_errors.BalancingError: Can't balance reaction by inv method"
```
The inverse algorithm can't handle this! Let's try others, like general pseudoinverse:
``` Python
>>> reaction.coefficients = reaction.balancer.gpinv()
[1954, 1854, 518, 1093, 1096, 55, 901, 898]
```
partial pseudoinverse:
``` Python
>>> reaction.coefficients = reaction.balancer.ppinv()
[39, 39, 13, 13, 11, 5, 16, 18]
```
and comb:
``` Python
>>> reaction.coefficients = reaction.balancer.comb()
[4, 5, 1, 1, 1, 1, 1, 3]
```

As we can see, we have got *four* different results (including 3 right ones) using *four* different algorithms. This is why the :meth:.ChemicalReaction.balancer.x() methods were implemented in the first place. We can, of course, get gpinv or ppinv data without intification:
``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
>>> reaction = ChemicalReaction(reaction_string, mode="balance", intify=False)

>>> reaction.coefficients = reaction.balancer.gpinv()
[1.4169688179840456, 1.3444525018129083, 0.37563451776649776, 0.7926033357505439, 0.7947788252356772, 0.03988397389412329, 0.6533720087019586, 0.6511965192168249]
```

!!! note
    Coefficients calculations by ChemicalReaction.balancer.x() methods works in all three modes (force, check and balance).

## Coefficients property setter
Unlike other properties of ChemicalFormula and ChemicalReaction, which are programmed as read-only cached properties without setters, the [coefficients][chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients] property can be set directly with appropriate list of coefficients:


``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="check")
>>> reaction.coefficients = [1, 5, 3, 3, 3, 3]

>>> reaction.coefficients
[1, 5, 3, 3, 3, 3]

>>> reaction.is_balanced
True
```

The [coefficients_validation()][chemsynthcalc.coefs.Coefficients.coefficients_validation] method checks setted coefficients **after setting them** in such properties as [normalized_coefficients][chemsynthcalc.chemical_reaction.ChemicalReaction.normalized_coefficients] and [masses][chemsynthcalc.chemical_reaction.ChemicalReaction.masses], therefore does not allow methods to calculate values with bad coefficients, for example:


``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
>>> reaction = ChemicalReaction(reaction_string, mode="check")
>>> reaction.coefficients = [1]

>>> reaction.coefficients
[1]

>>> reaction.is_balanced
False

>>> reaction.masses
chemsynthcalc.chem_errors.BadCoeffiecients: Number of coefficients should be equal to 6
```

This is one more way (along with direct input of coefficients in the reaction string for force and check modes) to input a custom set of coefficients for a particular reaction.

!!! note

    Coefficients setting works in all three modes (force, check and balance).

Now, when we covered *modes* and coefficients property calculations, we can get back to other arguments for reaction object creation.

## Target

The "target" parameter is the target chemical synthesis substance (i.e., whose mass is known in advance). The target choice is implemented as an integer pointer to formula index in products (right side of equation). There are, of course, other ways to do this (like explicitly input target as formula string), but integer index target was chose as method less prone to errors. Most of the time the target of the synthesis is the first product anyway (which is equal to 0 by default). We can set a target with object instantiation (the target is 1 which is FeO):

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
>>> reaction = ChemicalReaction(reaction_string, mode="balance", target=1)

>>> reaction.masses
[3.97359366, 0.28358172, 1.52731369, 1.0, 0.77944268, 0.12575572, 0.32138621, 0.5032771]
```

If we change the target, masses will obviously change too:

``` Python
>>> reaction = ChemicalReaction(reaction_string, mode="balance", target=2)

>>> reaction.masses
[5.09799345, 0.36382626, 1.95949455, 1.28296797, 1.0, 0.16134056, 0.41232821, 0.64568841]
```

The target can also be negative (limiting compoud - C in that case):

``` Python
>>> reaction = ChemicalReaction(reaction_string, mode="balance", target=-1)

>>> reaction.masses
[14.01216438, 1.0, 5.38579736, 3.52632041, 2.74856462, 0.443455, 1.13331074, 1.77471629]
```

## Target mass

The mass of target compound, the only mass that we know in advance before synthesis (in grams). We can change the target mass during the ChemicalReaction object instantiation (1 gram is the default):

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "Fe2O3+C=Fe3O4+FeO+Fe+Fe3C+CO+CO2"
>>> reaction = ChemicalReaction(reaction_string, mode="balance", target=1, target_mass=2)

>>> reaction.masses
[7.94718732, 0.56716344, 3.05462737, 2.0, 1.55888536, 0.25151145, 0.64277241, 1.0065542]
```

## ChemicalReaction properties

After the object initialization, we can access the ChemicalReaction properties:

* [reaction][chemsynthcalc.chemical_reaction.ChemicalReaction.reaction]
* [decomposed_reaction][chemsynthcalc.chemical_reaction.ChemicalReaction.decomposed_reaction]
* [_calculated_target][chemsynthcalc.chemical_reaction.ChemicalReaction._calculated_target]
* [chemformula_objs][chemsynthcalc.chemical_reaction.ChemicalReaction.chemformula_objs]
* [parsed_formulas][chemsynthcalc.chemical_reaction.ChemicalReaction.parsed_formulas]
* [matrix][chemsynthcalc.chemical_reaction.ChemicalReaction.matrix]
* [balancer][chemsynthcalc.chemical_reaction.ChemicalReaction.balancer]
* [molar_masses][chemsynthcalc.chemical_reaction.ChemicalReaction.molar_masses]
* [coefficients][chemsynthcalc.chemical_reaction.ChemicalReaction.coefficients]
* [normalized_coefficients][chemsynthcalc.chemical_reaction.ChemicalReaction.normalized_coefficients]
* [is_balanced][chemsynthcalc.chemical_reaction.ChemicalReaction.is_balanced]
* [final_reaction][chemsynthcalc.chemical_reaction.ChemicalReaction.final_reaction]
* [final_reaction_normalized][chemsynthcalc.chemical_reaction.ChemicalReaction.final_reaction_normalized]
* [masses][chemsynthcalc.chemical_reaction.ChemicalReaction.masses]
* [output_results][chemsynthcalc.chemical_reaction.ChemicalReaction.output_results]

## Output

A typical ChemicalReaction results output will look like this:

``` Python
>>> from chemsynthcalc import ChemicalReaction

>>> reaction_string = "KIO3+KI+H2SO4=I2+K2SO4+H2O"
>>> reaction = ChemicalReaction(reaction_string)
>>> reaction.print_results()

initial reaction: KIO3+KI+H2SO4=I2+K2SO4+H2O
reaction matrix:
 [[1. 1. 0. 0. 2. 0.]
 [1. 1. 0. 2. 0. 0.]
 [3. 0. 4. 0. 4. 1.]
 [0. 0. 2. 0. 0. 2.]
 [0. 0. 1. 0. 1. 0.]]
mode: balance
formulas: ['KIO3', 'KI', 'H2SO4', 'I2', 'K2SO4', 'H2O']
coefficients: [1, 5, 3, 3, 3, 3]
normalized coefficients: [0.33333333, 1.66666667, 1, 1, 1, 1]
algorithm: inverse
is balanced: True
final reaction: KIO3+5KI+3H2SO4=3I2+3K2SO4+3H2O
final reaction normalized: 0.33333333KIO3+1.66666667KI+H2SO4=I2+K2SO4+H2O
molar masses: [213.99947, 166.00247, 98.072, 253.80894, 174.252, 18.015]
target: I2
masses: [0.2811, 1.0901, 0.3864, 1.0, 0.6865, 0.071]
KIO3: M = 213.9995 g/mol, m = 0.2811 g
KI: M = 166.0025 g/mol, m = 1.0901 g
H2SO4: M = 98.0720 g/mol, m = 0.3864 g
I2: M = 253.8089 g/mol, m = 1.0000 g
K2SO4: M = 174.2520 g/mol, m = 0.6865 g
H2O: M = 18.0150 g/mol, m = 0.0710 g
```

One can output ChemicalReaction results using one of the 4 methods:

* [print_results][chemsynthcalc.chemical_reaction.ChemicalReaction.print_results]: print to stdout
* [to_txt][chemsynthcalc.chemical_reaction.ChemicalReaction.to_txt]: save as plain txt file
* [to_json][chemsynthcalc.chemical_reaction.ChemicalReaction.to_json]: serialization of output into an JSON object
* [to_json_file][chemsynthcalc.chemical_reaction.ChemicalReaction.to_json_file]: save as JSON file