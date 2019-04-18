RF Ablation example
===================

This example is based on the example described in the book "Modelling Organs, Tissues, Cells and Devices" by Socrates Dokos (https://www.springer.com/us/book/9783642548000).
Most of the parameters are taken from there. 
However, a full 3D model is presented here and a different geometry and electrode are used.

For the impedance boundary of the patch we assume for the real part of the admittance: conductivity @ 500 kHz divided by skin thickness (1.3mm).
For the imaginary part: 2 pi times 500 kHz times permittivity @ 500 kHz divided by skin thickness.
Values are taken from https://itis.swiss/virtual-population/tissue-properties/database/dielectric-properties/

The following studies are available:
* frequency-dependent part omitted: `parameters_RFA_direct.yml`
* stationary solution of the thermal problem: `parameters_RFAstationary.yml`
* frequency- and time-dependent study with Dirichlet BC on the patch: `parameters_RFA.yml` 
* frequency- and time-dependent study with Robin BC on the patch: `parameters_RFA_patch.yml` 
