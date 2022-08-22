# C-1-tricubic-treecode
A tricubic treecode method with global C^1 continuity

!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Authors:

  	Henry A. Boateng  (boateng@sfsu.edu) 
  	Department of Mathematics
  	San Francisco State University
  	San Francisco, CA
     
  	Svetlana Tlupova (tlupovs@farmingdale.edu)
  	Department of Mathematics
  	Farmingdale State College, SUNY
  	Farmingdale, NY
  
  This material is based upon work partially supported by the 
  National Science Foundation under Grant Nos. CHE-1800181 and DMS-2012371
  and U.S. Department of Energy, Office of Science, Office of Work- force 
  Development for Teachers and Scientists (WDTS) under the Visiting Faculty Program.
  
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


   NOTE: Please include the following references in any work that
         utilizes this code:
		 
        (1) Boateng. H. A., Tlupova, S.: A treecode algorithm based on 
            tricubic interpolation
            submitted, (2022)  
		 

Summary of files :
------------------

      rand*.txt       : Files containing the coordinates of the particles
                        distributed in [-5.0,5.0]x[-5.0,5.0]x[0.0,10.0].
                        The files rand_10000.txt, rand_80000.txt and rand_640000.txt
                        correspond to 10000, 80000, and 640000 particles respectively.
                        
      lambda*.txt      : Files containing the weights of the particles in the
                         corresponding rand*.txt
                        
                      

!!!!!!!
! .cpp and .h FILES 
      direct_sum.cpp : C++ file for exact computation of the potential and field. 
      
      Tricubic.cpp   : Main C++ program for the tricubic treecode.
      
                       Both direct_sum.cpp and Tricubic.cpp depend on the 
                       following helper files (utilities.h, utilities.cpp, kernel_utils.h,
                       kernel_utils.cpp, tree_utils.h, tree_utils.cpp, tricubic_utils.h,
                       tricubic_utils.cpp)
      
      makefile       : A makefile for compiling the direct sum and tricubic treecode. 
                       It produces the executables direct_sum and Tricubic

      input_params.txt   : Input file for the both direct_sum and Tricubic


Input for direct_sum and Tricubic :
-----------------------------------

      input_params.txt specifies the following required options in
      the given order:
      
#// Treecode method to be used: Particle-Cluster, Cluster-Particle

	Particle-Cluster    #// (This option picks particle-cluster)
 
#// Kernel name below: Coulomb, ScreenedCoulomb, RSEwaldSum

	RSEwaldSum // (This option picks the  real space Ewald Sum kernel)
 
// number of particles

	10000            
 
// N0, leaf size 

	2000     // Pick N0 = 2000
 
// parameter to choose the type of MAC (0 - Regular MAC, 1 - Spherical MAC)

	0   // Use Regular MAC
 
// theta, MAC parameter (Should be less than 1 for Regular MAC, ideally greater than 1 for Spherical MAC)

	0.5  //theta = 0.5
 
// 1 if finite differences should be used, 0 otherwise

	0  // Use analytical derivatives
 
// mgrid, any value when finite differences are not used. Do not delete.

	20
 
// hgrid, any value when finite differences are not used. Do not delete.

	0.01

Output for the executable direct_sum :
-------------------------------------

The output is the file exact_sum_kernelName_Nnumberofparticles. 
Example: The output for the Coulomb sum with 10000 particles is exact_sum_Coulomb_N10000

For a system of size N, the file has N+1 lines. The first N lines correspond to data for
the N particles. Each of the first N lines has 5 entries:

	Particle number, particle potential, x-component of gradient of potential, y-component of gradient of potential, z-component of gradient of potential
   
The N+1 line is the time the computation took in 10^(-2) seconds (i.e. centiseconds)


Output for the executable Tricubic  :
-----------------------------------

Running Tricubic generates two files : output.txt and tricubic_sum_kernalName_LM_Nnumberofparticles

The file output.txt has the following entries:

	Number of particles (N_cube), Maximum number of particles in a leaf (N0), MAC (theta), relative 2-norm error in potential (E), relative 2-norm error in field (FE), relative 1-norm error in potential (eL1), relative 1-norm error in field, time for direct sum in centiseconds, treecode time in centiseconds

For the file tricubic_sum_kernalName_LM_Nnumberofparticles:
    The first line has the data:
    cpu time for direct sum (in centiseconds), cpu time for treecode (in centiseconds)
    
    Following the first line are direct sum and the treecode approximation for the potential and fields in the sample format:
    
    particle number,  exact potential, treecode approximation of potential
    x-component of exact gradient of potential, y-component of exact gradient of potential, z-component of exact gradient of potential
    tree approximation of x-component of gradient of potential, treecode approximation of y-component of gradient of potential, treecode approximation of of z-component of gradient of potential


