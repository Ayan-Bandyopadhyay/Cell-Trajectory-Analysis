Cell-Trajectory-Analysis
========================

##Synopsis
**Cell-Trajectory-Analysis** is a set of functions for cell tracking and biological trajectory analysis, to be used on video
taken with fluorescence microscopy. This functions can be used to threshold video to detect cells marked with green fluorescence protein, generate trajectories for each cell, quantify each trajectory, and classify each trajectory using the k-nearest neighbors statistical learning model. This package is dependent on two other packages, namely flowcatchR and traj.

##Documentation 

###fitModelKNN_CV
**Description**
####Usage

####Arguments
####Value
####Author
####Examples



###generateTraj
####Description
####Usage
generateTraj(particles, L = 26, R = 3, epsilon1 = 0, epsilon2 = 0,
  lambda1 = 1, lambda2 = 0, penaltyFunction = penaltyFunctionGenerator(),
  include.area = FALSE, frames)
####Arguments
particles	          A ParticleSet object
L	                  Maximum number of pixels an object can move in two consecutive frames
R	                  Linkrange, i.e. the number of consecutive frames to search for potential candidate links
epsilon1	          A numeric value, to be used in the formula. Jitter for allowing angular displacements
epsilon2	          A numeric value, to be used in the formula. Jitter for allowing spatial displacements
lambda1	            A numeric value. Multiplicative factor for the penalty function
lambda2	            A numeric value. Multiplicative factor applied to the angular displacement
penaltyFunction     A function structured in such a way to be applied as penalty function in the linking
include.area	      Logical, whether to include also area change of the particles in the cost function calculation
frames	            The Frames object that the ParticleSet object is derived from
verboseOutput	      Logical, whether the output should report additional intermediate steps. For debugging use mainly.
prog	              Logical, whether the a progress bar should be shown during the tracking phase
include.intensity	  Logical, whether to include also intensity change of the particles in the cost function calculation
####Value
A TrajectorySet object
####Author
Ayan Bandyopadhyay, Bellarmine College Prep
####Examples
```r
library(flowcatchR)
platelets <-particles(channel.Frames(MesenteriumSubset,"red"))
trajSet <- generateTraj(platelets,
                       L=26, R=3,
                       epsilon1=0, epsilon2=0,
                       lambda1=1, lambda2=0,
                       penaltyFunction=penaltyFunctionGenerator(),
                       include.area=FALSE, MesenteriumSubset)
```

###greenThresh
####Description
####Usage
####Arguments
####Value
####Author
####Examples


###trajMeasures
####Description
####Usage
####Arguments
####Value
####Author
####Examples
