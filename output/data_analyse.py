import numpy as np
import matplotlib.pyplot as plt
import os
print(os.getcwd())
potential_energy, kinetic_energy, total_energy, ps = np.loadtxt("outfile.txt", unpack=True)

plt.figure()
plt.plot(ps, total_energy, '-')
plt.xlabel('ps')
plt.ylabel('total energy')
plt.savefig('total_energy.png')

plt.figure()
plt.plot(ps, kinetic_energy, '-')
plt.xlabel('ps')
plt.ylabel('kinetic energy')
plt.savefig('kinetic_energy.png')

plt.figure()
plt.plot(ps, potential_energy, '-')
plt.xlabel('ps')
plt.ylabel('potential energy')
plt.savefig('potential_energy.png')


