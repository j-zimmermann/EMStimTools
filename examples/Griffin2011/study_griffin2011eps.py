from emstimtools import Simulation

model = Simulation('parameters_griffin2011eps.yml')
model.set_log_level('DEBUG')
model.run()
print("success")
