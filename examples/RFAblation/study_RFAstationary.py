from emstimtools import Simulation

model = Simulation('parameters_RFAstationary.yml')
model.set_log_level('DEBUG')
model.run()
print("success")
