from emstimtools import Simulation

model = Simulation('parameters_mobini2017.yml')
model.set_log_level('DEBUG')
model.run()
print("success")

