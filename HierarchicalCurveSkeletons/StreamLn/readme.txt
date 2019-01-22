StreamLn module
---------------
	Given a force field, a set of critical points, 
and other seeds (high divergence, high curvature), finds a set of skeleton 
points by following the streamlines originating at critical points and the 
seeds.

** NOT USED - StreamLn-Search.cpp (includes detection of critical points)
   ------------------------------
	Finds the critical points by interpolating the force field on a 
10x10x10 grid in each voxel and looking for points where the force is smaller
than a certain threshold.
	Main Problem: - how to choose the threshold so that it works for any 
dataset.
	Follows streamlines originating at all detected critical points, in 
the direction pointed by the force vector at that point. If the force vector 
is exactly 0 ... ? ;)
	Uses a flag array of the same size as the volume. A flag is set to 1 
if a skeleton point lies within that voxel.

** NOT USED - StreamLn-flags.cpp
   -----------------------------
   	Based on StreamLn-Search.cpp.
   	Does not include detection of critical points. Critical points are 
input to the module.
	Follows streamlines starting only at saddle points and stops if the 
current voxel is an already visited voxel (the flag is set), or if a neighbor 
of the current voxel is already visited.

** NOT USED - StreamLn-allInOne.cpp
   --------------------------------
	Critical points are input to the module.
   	Doesn't use a flags array.
   	Follows streamlines starting at:
   		- critical points of type saddle (in the direction of the 
   			positive eigenvectors of the Jacobian evaluated at the 
   			critical point, and the exact opposite direction)
and goes on until close to a point that is already part of the skeleton. Close
means the Manhattan distance between the 2 points is smaller than a certain 
threshold.
	The whole skeleton is computed in one function (connects critical 
points, and high div points).


** CURRENT VERSION - StreamLn.cpp 
   ------------------------------
   	Based on StreamLn-allInOne.cpp. Separates the code into a set of 
functions. One of them computes the basic skeleton, connecting only the 
critical points. Another one, adds the streamlines seeded at high divergence 
points. Another one adds the streamlines starting at boundary seeds.
