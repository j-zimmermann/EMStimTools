from emstimtools import Simulation

model = Simulation('parametersRobinBCNaturalConvection.yml')
model.set_log_level('DEBUG')
model.run()
print("success")
