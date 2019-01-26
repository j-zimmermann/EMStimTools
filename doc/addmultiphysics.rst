How to add a new multiphysics interface
---------------------------------------

#. Test the study with a pure FEniCS script.
#. When the study is working for an example and validated, separate it into different part:
    #. Boundary Conditions: check if the currently implemented interface is sufficient. If so, generalise your model w.r.t to the BCs.
    #. Implement material parameters: take them and use the :meth:`emstimtools.fenics.fenics.Fenics._set_material_constant` to set them properly. Pay attention to implementing the assignment properly! The vector containing the values is not assigned automatically by just taking the name of the parameter since :func:`setattr` does not work well with nested objects.   
    #. take care of file output and derived quantities.
