Example for Capacitive Coupling
===============================

The set-up is taken from a paper by Griffin et al. (dx.doi.org/10.1371/journal.pone.0023404).

Run
---

Use 

.. code::

	python3 study_griffin2011.py

or another study. HF refers to 60kHz, eps to eps_r=3000.

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
Go for instead with a faster solution in paraview (.pvsm file in paraview directory).
