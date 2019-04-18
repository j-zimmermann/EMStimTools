from emstimtools import Simulation

model = Simulation('parametersRobinBCNaturalConvectionTimeDependent.yml')
model.set_log_level('DEBUG')
model.run()
print("success")
