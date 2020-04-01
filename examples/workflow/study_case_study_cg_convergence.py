import yaml
import numpy as np
import matplotlib.pyplot as plt
from emstimtools import Simulation


def prepare_input_dict(max_size, min_size, data):
    data['salome_parameters']['max_size'] = max_size
    data['salome_parameters']['min_size'] = min_size
    data['mesh'] = {}
    data['mesh']['transformation'] = {}
    data['mesh']['transformation']['scale'] = 1e3
    with open('parameters_case_study_cg_convergence.yml', 'w') as stream:
        yaml.dump(data, stream)
 
max_size = 0.024
min_size = 0.001

iterations = 7 
max_list = []
min_list = []
dof_list = []
current1_list = []
current2_list = []
field_list = []
error_list = [np.nan]

with open('parameters_case_study_cg_convergence_master.yml', 'r') as stream:
    input_dict = yaml.safe_load(stream)
input_dict_file = 'parameters_case_study_cg_convergence.yml'
for i in range(iterations):
    prepare_input_dict(max_size, min_size, input_dict)
    model = Simulation(input_dict_file)
    model.run()
    if i > 0:
        error = model.fenics_study.get_solution_error(old_solution, save_difference=True)
        l2error = model.fenics_study.get_solution_l2error(old_solution)
        error_list.append((error, l2error))

    dofs = model.fenics_study.V.dim()
    model.fenics_study.get_field()
    current1 = model.fenics_study.compute_current("Contact1")
    current2 = model.fenics_study.compute_current("Contact1")
    old_solution = model.fenics_study.u.copy(deepcopy=True)
    dof_list.append(dofs)
    max_list.append(max_size)
    min_list.append(min_size)
    current1_list.append(current1)
    current2_list.append(current2)
    if i > 0:
        print("Current error: ", abs(current1_list[-1] - current1_list[-2]), flush=True)
    field_list.append(model.fenics_study.Efield(0.0, 0.0, 0.001))
    max_size = max_size / 2.
    min_size = min_size / 2.
    print("##############", flush=True)
    print("Current state:", flush=True)
    print("##############", flush=True)
    print(dof_list, max_list, min_list, current1_list, current2_list, field_list, error_list, flush=True)
    print("##############", flush=True)
    print("##############", flush=True)

field_norms = [np.linalg.norm(np.array(f)) for f in field_list]

plt.xscale('log')
plt.xlabel("DOFs")
plt.ylabel("Current [mA]")
plt.plot(dof_list, current1_list)
plt.savefig("current.pdf")
plt.close()

plt.xscale('log')
plt.xlabel("DOFs")
plt.ylabel("Electric field [V/m]")
plt.plot(dof_list, 1e3 * np.array(field_norms))
plt.savefig("field.pdf")
plt.close()
print ("done")
