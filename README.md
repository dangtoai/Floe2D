# Floe2D
We studie the fracture by collision at the scale of ice floe size.
We suppose that each floe is an elastic, homogenous, isotropic material so that it is legitime to use a mass-spring network to compute its deformation after shock.

We focus on the deformation field applied on a floe after a certain time so that we can use its trace to compute fracture thanks to the theory of Francfort and Marigo. 


Install necessary modules: pip3 install -r requirements.txt

## how to use:

to obtain a masses-springs network of 1 geometry 
```
python3 generate_network.py 
```
then enter: 3 

to obtain a simulation of a collision on the network
```
python3 deformation_computation.py 
```
enter: 20 20 
enter: -30 20

to generate a mesh from GMSH to solve the griffith energy
```
python3 generated_geofile.py 
```
then enter: 3

to run a fracture computation on the geometry above: 
```
python3 BoundarySolver.py 
```


