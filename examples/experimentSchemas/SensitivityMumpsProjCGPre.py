import pandas as pd
import numpy as np
import itertools
import yaml
from emstimtools import Simulation


# input factors and factor levels
# TODO insert parameter names and base and sensitivity values here
dat = [('conductivity_medium', 1.0, 1.5),
       ('permittivity_medium', 60, 80),
       ('permittivity_dish', 2.5, 4.0)]
inputs_df = pd.DataFrame(dat, columns=['factor', 'baseCase', 'sensitivityCase'])
inputs_df = inputs_df.set_index(['factor'])
print(inputs_df)

# get number of factors
k = len(dat)

# 2 level full factorial experiment design, -1 base case, 1 sensitivity case
doe_df = pd.DataFrame(list(itertools.product([-1, 1], repeat=k)))
print(doe_df)

# generate configurations based on experiment design
configs_df = doe_df.copy(deep=True)
for i in range(0, 2**k):
    for j in range(0, k):
        if doe_df.loc[i, j] < 0:
            # set base value
            configs_df.iloc[i, j] = inputs_df.iloc[j, 0]
        else:
            # set sensitivity value
            configs_df.iloc[i, j] = inputs_df.iloc[j, 1]
print(configs_df)

# initialise parameters
with open('parameters_griffin2011meshmaster.yml', 'r') as stream:
    data = yaml.load(stream)


# run the model for all configurations and get results
results = []

for i in range(0, 2**k):
    # update parameters
    data['conductivity']['medium'] = float(configs_df[0][i])
    data['permittivity']['medium'] = float(configs_df[1][i])
    data['permittivity']['dish'] = float(configs_df[2][i])
    data['solver']['linear_solver'] = 'mumps'
    data['properties']['project_solver'] = 'cg'
    data['properties']['project_preconditioner'] = 'petsc_amg'

    with open('parameters_griffin2011mesh.yml', 'w') as stream:
        yaml.dump(data, stream)

    model = Simulation('parameters_griffin2011mesh.yml')
    model.run()
    results.append(model.fenics_study.normEr(0.05, 0.05, 0.01101))

results_df = pd.DataFrame(results)
print(results_df)


# Compute the individual effect of a factor on the response
individual_effects = list(itertools.repeat(0, k))
total_effects = list(itertools.repeat(0, k))
interaction_effects = list(itertools.repeat(0, k))
# Go through configurations
for row in range(0, 2 ** k):
    indList = []
    # Get the change in the output for this configuration
    change = results_df.iloc[row, 0] - results_df.iloc[0, 0]
    # Check the factor values
    for fac in range(0, k):
        # if sensitivity case
        if doe_df.iloc[row, fac] > 0:
            # add change to the total effects of the factor f
            total_effects[fac] = total_effects[fac] + change
            # store in helper list
            indList.append(fac)
    if len(indList) == 1:
        # add as individual effect
        individual_effects[indList[0]] = change

# interaction effects
interaction_effects = np.subtract(total_effects, individual_effects)


print('Individual Effects ', individual_effects)

print('Total Effects ', total_effects)

print('Interaction Effects ', interaction_effects)
