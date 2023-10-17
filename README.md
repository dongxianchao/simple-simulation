# simple-simulation
This is a simple program I wrote for molecular dynamics simulations, intended solely for learning purposes.
Now it supports calculations using the LJ potential, EAM potential and EAM/alloy potential. For specific input formats, you can refer to my 'incar' file.

Now I have updated the NVT ensemble and I inserted 2 method to simulate the temperature. You can use "ber" or "bdp" to set the mode, ber refers to Berendsen thermostat and bdp refers to Bussi-Donadio-Parrinello thermostat. 


TO DO:
I'm implementing a thermostat and barostat algorithm to control pressure during the simulation.
