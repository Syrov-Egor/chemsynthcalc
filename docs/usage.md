Let's say that we need to prepare 3 grams of [YBCO](https://en.wikipedia.org/wiki/Yttrium_barium_copper_oxide) by solid-state synthesis from respective carbonates.

## Reaction string

``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"
```

## ChemicalReaction object

Now, we can create a chemical reaction object of [ChemicalReaction][chemsynthcalc.chemical_reaction.ChemicalReaction] class, which will be used in the calculation. We need to specify arguments for our particular case:
``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"

reaction = ChemicalReaction(
    reaction = reaction_string, # our reaction string
    target = 0, # number of target compound in product list
    target_mass = 3, # desired mass of target compound,
    mode = "balance" # mode of coefficients calculations,
)
```

## Calculation and output

Now, to perform automatic calculation, all we need to do is to put [print_results()][chemsynthcalc.chemical_reaction.ChemicalReaction.print_results] method:
``` Python
from chemsynthcalc import ChemicalReaction

reaction_string = "BaCO3 + Y2(CO3)3 + CuCO3 + O2 → YBa2Cu3O7 + CO2"

reaction = ChemicalReaction(
    reaction = reaction_string, # our reaction string
    target = 0, # number of target compound in product list
    target_mass = 3, # desired mass of target compound,
    mode = "balance" # mode of coefficients calculations,
)

reaction.print_results(print_rounding_order=4)
# assuming we use analytical balances with 4 digits presicion
```

And we get our ouput in the terminal:

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

Thus, we got all masses ready for our planned synthesis!