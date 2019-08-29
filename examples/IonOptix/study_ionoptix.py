from emstimtools import Simulation

model = Simulation('parameters_ionoptix.yml')
model.set_log_level('DEBUG')
model.run()
print("success")

