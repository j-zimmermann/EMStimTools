from emstimtools import Simulation

# dielectric heating
# model = Simulation('parameters_RFA.yml')
# model.set_log_level('DEBUG')
# model.run()

# direct Joule heating (no frequency dependence)
model = Simulation('parameters_RFA_direct.yml')
model.set_log_level('DEBUG')
model.run()

# dielectric heating with surface impedance on patch
# model = Simulation('parameters_RFA_patch.yml')
# model.set_log_level('DEBUG')
# model.run()

print("success")
