from emstimtools import Simulation

model = Simulation('parameters_case_study_cg_mesh.yml')
model.set_log_level('DEBUG')
model.run()
print("success")
