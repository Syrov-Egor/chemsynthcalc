#Installing

##Python 3.10 or later

chemsynthacalc requires [Python 3.10](https://www.python.org/downloads/) or later.

##Installation

Install from [pypi](https://pypi.org/):

``` bash
pip install chemsynthcalc
```


##NumPy and SciPy

NumPy and SciPy are requirements for fast matrix operations for reaction balancing. They will be installed automatically by pip if they are not already installed.

##But why chemsynthcalc is using SciPy?

In short: to ensure consistent results across all platforms. See [Note on SciPy][chemsynthcalc.balancing_algos.BalancingAlgorithms].