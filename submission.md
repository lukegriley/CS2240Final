## FEM (final submission)


### Output
See video:
https://drive.google.com/file/d/1hkqQdZKfSBudEc16rXedYG5Sh09-pvr1/preview

## Extra Features
- RK4 implicit integration. A total of 4 euler deriv evals are taken per timestep, leading to more stable simulation.
- Parallelization (shown in video): I used OpenMP in three places:
  - Looping through each node to detect collisions.
  - Looping through each node to zero forces + apply gravity.
  - Looping through each tet to apply internal forces (this saved the most time).
- Interactivity:
  - Push (by pressing Z): Applies a force to every node in the mesh in the opposite direction of the camera. The force is proportional to the dot product of the camera's look and camera's direction to the node position, minus cos(pi/4).
  - Pull (by pressing X): Applies a force to every node in the mesh in the direction of the camera. The force on each node has the same proportions as in push.
- Frame-by-frame stepping (by pressing P): runs one update cycle, when the simulation is paused, to allow for precise timestepping.
- Octocat mesh: I could only finding .msh tet files online, consisting of binary data. I imported one into the open source Gmsh software, which could export to .mesh. Although it was no longer binary data, it was still not in the correct format for the mesh loader header. I created a python script (see meshconverter.py) to convert into the acceptable format. From there, I had a workflow to import any .msh file, such as the octopus-cat aka OCTOCAT mesh which I downloaded from the fTetWild repo.

- Unfinished: plasticity using Green's strain.


### Design Choices
```Initialization```
Most of the FEM logic was done in system.cpp, which is initialized by Simulator simple vertex position data and tet-vertex indices. 
It creates Tet and Node structs from this. Important information is derived from the undeformed shape upon initialization. I use the material space coordinates to create the beta matrix, which is later used to derive the deformation gradient.
I also, for the only time during the simulation, calculate the areas and normals of each tet's face, as well as node masses, at the start.

```Update Cycle```
I implemented the classic Runge-Kutta method in my System's update() function. I start by copying each node's position and velocity into a single oldPosVel vector. 
Each derivative is added to this vector. The subsequent derivative is evaluated with node positions as the old positions added to the previous derivative. 

Inside the deriv eval function, I zero out all node forces and apply gravity. If the user had press the force push button, I had the push/pull force as well. After, I call applyInternalForces(), which calculates the deformation gradient, strain, stress, and viscous damping. 

```Collisions```
Only after the implicit integration is complete and the final positions are updated, I calculate and apply collision "forces". Here, I project the colliding node onto the collider's surface, decompose its velocity into a normal and tangential component, apply friction, and invert the normal component. 
I calculate collisions after RK4 so the velocity is not being inverted and reduced multiple times per update cycle.



### Collaboration/References
Lecture slides

RK4: https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

OpenMP: http://supercomputingblog.com/openmp/tutorial-parallel-for-loops-with-openmp/

https://en.wikipedia.org/wiki/Tangential_and_normal_components

Meshes: https://github.com/wildmeshing/fTetWild/tree/master

Plasticity: http://graphics.berkeley.edu/papers/Obrien-GMA-2002-08/Obrien-GMA-2002-08.pdf

Vector hash: https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/


### Known Bugs
Meshes that aren't smartly constructed are very unstable, as can be seen with the octocat mesh in the video.  
