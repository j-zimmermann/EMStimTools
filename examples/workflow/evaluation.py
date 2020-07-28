import numpy as np
import matplotlib.pyplot as plt

ndofs = [10845, 38027, 132482, 141123, 212993, 756790, 5125719]
max_size = [0.024, 0.012, 0.006, 0.003, 0.0015, 0.00075, 0.000375]
min_size = [0.001, 0.0005, 0.00025, 0.000125, 6.25e-05, 3.125e-05, 1.5625e-05] 
current1 = [-7.249135556530399, -10.285699768588968, -11.642418005757829, -11.661499703428715, -11.655092793422824, -11.64002780441075, -11.670704461364625] 
current2 = [-7.249135556530399, -10.285699768588968, -11.642418005757829, -11.661499703428715, -11.655092793422824, -11.64002780441075, -11.670704461364625] 
field = [np.array([ 1.40638964e-05,  3.22568829e-02, -1.40824838e-04]), np.array([3.58654541e-05, 3.32119444e-02, 1.71986016e-05]), np.array([ 6.56293087e-06,  3.30352692e-02, -2.64033628e-06]), np.array([-3.35717911e-06,  3.30305701e-02,  9.55175351e-05]), np.array([-2.64203286e-06,  3.28287092e-02,  4.64169762e-05]), np.array([5.73693438e-07, 3.28292303e-02, 9.91321045e-06]), np.array([ 1.31449348e-08,  3.28057970e-02, -1.63257845e-06])] 
errors = [(np.nan, np.nan), (1.6352517938125766, 0.5554247915574236), (1.0936455395460682, 0.1738828809554651), (0.10524639054394395, 0.04462514356457091), (0.07226088787267955, 0.02156715354139777), (0.09137641329711556, 0.012712789533963504), (0.19327818729008178, 0.013982047525451142)]


errorsl2 = [x[1] for x in errors]

current_error = [np.nan] 
for i in range(1, len(current1)):
    current_error.append(np.abs(current1[i-1] - current1[i]))

field_norms = [1e3 * np.linalg.norm(f) for f in field]
ticks = ["({:.2e}, {:.2e})".format(min_size[i], max_size[i]) for i in range(len(max_size))] 

plt.xscale('log')
plt.xlabel("DOFs")
plt.ylabel("current [mA]")
plt.plot(ndofs, np.array(current1), '+-')
plt.savefig("current1_conv.pdf", dpi=300)
plt.show()

plt.xscale('log')
plt.yscale('log')
plt.xlabel("DOFs")
plt.ylabel("current error [mA]")
plt.plot(ndofs, np.array(current_error), '+-')
plt.savefig("current1error_conv.pdf", dpi=300)
plt.show()

plt.yscale('log')
plt.xlabel("SALOME configuration (min_size, max_size)")
plt.xticks(np.arange(len(max_size)), ticks, rotation=40)
plt.ylabel("current error [mA]")
plt.plot(np.array(current_error), '+-')
plt.tight_layout()
plt.savefig("current1error_conv_SALOME.pdf", dpi=300)
plt.show()

plt.xscale('log')
plt.xlabel("DOFs")
plt.ylabel("electric field strength [V/m]")
plt.plot(ndofs, field_norms, '+-')
plt.savefig("field_conv.pdf", dpi=300)
plt.show()

plt.xlabel("SALOME configuration (min_size, max_size)")
plt.ylabel("electric field strength [V/m]")
plt.xticks(np.arange(len(max_size)), ticks, rotation=40)
plt.plot(field_norms, '+-')
plt.tight_layout()
plt.savefig("field_conv_SALOME.pdf", dpi=300)
plt.show()

plt.xlabel("SALOME configuration (min_size, max_size)")
plt.ylabel("current [mA]")
plt.xticks(np.arange(len(max_size)), ticks, rotation=40)
plt.plot(np.array(current1), '+-')
plt.tight_layout()
plt.savefig("current1_conv_SALOME.pdf", dpi=300)

plt.show()

plt.xscale('log')
plt.yscale('log')
plt.xlabel("DOFs")
plt.ylabel("potential [V], L2 error")
plt.plot(ndofs, errorsl2, '+-')
plt.tight_layout()
plt.savefig("errors.pdf", dpi=300)
plt.show()

plt.yscale('log')
plt.xlabel("SALOME configuration (min_size, max_size)")
plt.ylabel("potential [V], L2 error")
plt.xticks(np.arange(len(max_size)), ticks, rotation=40)
plt.plot(errorsl2, '+-')
plt.tight_layout()
plt.savefig("errors_SALOME.pdf", dpi=300)
plt.show()

#ndofs.pop(-1)
#errorsl2.pop(-1)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("DOFs")
plt.ylabel("potential [V], L2 error")
plt.plot(ndofs, errorsl2, '+-')
plt.tight_layout()
plt.savefig("errors_corr.pdf", dpi=300)
plt.show()
print("Error rate:\n")
for i in range(2, len(ndofs)):
    print(np.log(errorsl2[i] / errorsl2[i-1]) / np.log(ndofs[i] / ndofs[i-1]))
