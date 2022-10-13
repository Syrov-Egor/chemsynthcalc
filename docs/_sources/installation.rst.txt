Installing
============

Python 3.8 or later
-------------------

chemsynthacalc requires `Python 3.8 <https://www.python.org/downloads/>`_ or later.

Installation
------------

Install from `pypi <https://pypi.org/project/chemsynthcalc/>`_::

  pip install chemsynthcalc

After installation, run test to make sure everything works properly::

  import chemsynthcalc

  chemsynthcalc.run_test()

NumPy and SciPy
---------------

`NumPy <https://numpy.org/>`_ and `SciPy <https://scipy.org/>`_ are requirements 
for fast matrix operations for reaction balancing. They will be installed automatically
by pip if they are not already installed.

But why chemsynthcalc is using SciPy?
--------------------------------------

In short: to ensure consistent results across all platforms.
See :class:`chemsynthcalc.reaction_balance.Balancer` note for the full answer.