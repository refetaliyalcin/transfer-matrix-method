# transfer-matrix-method

Solver belongs to Shawn Divitt from ETH Zurich. Before finding this code, I was using https://github.com/sbyrnes321/tmm that one which written in Python. I couldn't find a robust solver for MATLAB which can handle different incident angles and wavelengths, then I find this robust code. I only made a small change which I implemented from https://github.com/sbyrnes321/tmm to improve stability for high absorptive layers.

Solves reflectance, transmittance and absorptance of plane parallel multilayered coatings at different angles.

Solver is jreftran_rt.m however main code to run is main.m which loops the solver over wavelengths and incident angles.
