import yaml
import numpy as np
import matplotlib.pyplot as plt
from emstimtools import Simulation

initial_size = 0.012
min_size = 0.0005

size = initial_size
iterations = 7 
max_list = []
min_list = []
dof_list = []
current1_list = []
current2_list = []
field_list = []

with open('parameters_case_study_cg_convergence_master.yml', 'r') as stream:
    data = yaml.load(stream)

for i in range(iterations):
    data['salome_parameters']['max_size'] = size
    data['salome_parameters']['min_size'] = min_size
    with open('parameters_case_study_cg_convergence.yml', 'w') as stream:
        yaml.dump(data, stream)
    model = Simulation('parameters_case_study_cg_convergence.yml')
    # model.set_log_level('DEBUG')
    model.run()
    dofs = model.fenics_study.V.dim()
    model.fenics_study.get_field()
    current1 = model.fenics_study.compute_current("Contact1")
    current2 = model.fenics_study.compute_current("Contact1")
    dof_list.append(dofs)
    max_list.append(size)
    min_list.append(min_size)
    current1_list.append(current1)
    current2_list.append(current2)
    field_list.append(model.fenics_study.Efield(0.0, 0.0, 0.001))
    size = size / 2.
    min_size = min_size / 2.
    print("##############", flush=True)
    print("Current state:", flush=True)
    print("##############", flush=True)
    print(dof_list, max_list, min_list, current1_list, current2_list, field_list, flush=True)
    print("##############", flush=True)
    print("##############", flush=True)

field_norms = [np.linalg.norm(np.array(f)) for f in field_list]

plt.xscale('log')
plt.plot(dof_list, current1_list)
plt.savefig("current.pdf")
plt.close()

plt.xscale('log')
plt.plot(dof_list, field_norms)
plt.savefig("field.pdf")
plt.close()

