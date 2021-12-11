how to run:
open hydro.m
uncomment line 13 or 14 depending if you are testing the Newtonian or Lagrangian dynamics
run hydro.m

to change parameters, either directly change the parameters.m file, 
or pass in a modified parameter struct to hydro.m as a function input

the lagrangian symbolic math is done in lagrangian.mlx and is exported to generated_accels.m