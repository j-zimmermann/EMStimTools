import uncertainpy as un
import chaospy as cp
import yaml
from emstimtools import Simulation


def griffinmodel(conductivity_medium, permittivity_medium, permittivity_dish):

    with open('parameters_griffin2011meshmaster.yml', 'r') as stream:
        data = yaml.load(stream)

    # update parameters
    data['conductivity']['medium'] = float(conductivity_medium)
    data['permittivity']['medium'] = float(permittivity_medium)
    data['permittivity']['dish'] = float(permittivity_dish)
    data['solver']['linear_solver'] = 'mumps'

    with open('parameters_griffin2011mesh.yml', 'w') as stream:
        yaml.dump(data, stream)

    model = Simulation('parameters_griffin2011mesh.yml')
    model.set_log_level("DEBUG")
    model.run()

    field = model.fenics_study.normEr(0.05, 0.05, 0.01101)
    print("Solved for conductivity of medium: ", data['conductivity']['medium'])
    print("Solved for permittivity of medium: ", data['permittivity']['medium'])
    print("Solved for conductivity of dish: ", data['conductivity']['dish'])
    print("Got field: ", field)

    return 1, field


parameters = {'conductivity_medium': cp.Uniform(1.0, 1.5),
              'permittivity_medium': cp.Uniform(60, 80),
              'permittivity_dish': cp.Uniform(2.5, 4.0)
              }

uq_model = un.Model(griffinmodel, labels=[r"Electric Field Strength [V/m]"])
UQ = un.UncertaintyQuantification(model=uq_model, parameters=parameters, CPUs=None)

data = UQ.quantify(seed=42, method="pc")
