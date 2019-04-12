Example for Capacitive Coupling and Sensitivity Analysis
========================================================

The set-up is taken from a paper by Griffin et al. (dx.doi.org/10.1371/journal.pone.0023404).
This study is a first example iof how a sensitivity analysis could look like.

Run
---

Use 

.. code::

	python3 study_griffin2011sensitivity.py


Postprocessing
--------------

If you don't have much time, remove the following lines from the parameter YAML-file:

.. code::

	properties:
        	E-Field : yes
        	E-Field-norm : yes
        	project_element : DG
        	project_degree : 1

This is one the hand a projection that yields a correct field and norm value.
However, it takes very long.
Instead, you could use paraview filters.
