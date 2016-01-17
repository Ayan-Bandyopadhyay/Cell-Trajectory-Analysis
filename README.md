Cell-Trajectory-Analysis
========================

##Synopsis
**Cell-Trajectory-Analysis** is a set of functions for cell tracking and biological trajectory analysis, to be used on video
taken with fluorescence microscopy. This functions can be used to threshold video to detect cells marked with green fluorescence protein, generate trajectories for each cell, quantify each trajectory, and classify each trajectory using the k-nearest neighbors statistical learning model. This package is dependent on two other packages, namely flowcatchR and traj.

##Documentation 

###fitModelKNN_CV
**Description**
#
Fits trajectory data into a K-nearest neighbors classifier using leave one out cross validation.
#
**Usage**
#
fitModelKNN_CV(trajDataFrame, labelVector, kVal)
#
**Arguments**
#
trajDataFrame &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A data frame with measures for each trajectory. This can be produced by function trajMeasures
#
labelVector	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A vector of labels for each trajectory. Its length must be equal to the number of rows in trajDataFrame
#
kVal &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;	The number of neighbors used for classification
#
**Value**
#
A vector of classified labels for each trajectory
#
**Author**
#
Ayan Bandyopadhyay, Bellarmine College Prep
#
**Examples**
```r
data <- as.data.frame(matrix(1:4,nrow = 2,ncol = 2))
labels <- c("live","dead")
classifierKNN_CV<- fitModelKNN_CV(data,labels,3)
```

###generateTraj
#
**Description**
#
Returns a TrajectorySet so that the y value is the distance from the bottom of the Frames object, not from the top.
#
**Usage**
#
generateTraj(particles, L = 26, R = 3, epsilon1 = 0, epsilon2 = 0,
  lambda1 = 1, lambda2 = 0, penaltyFunction = penaltyFunctionGenerator(),
  include.area = FALSE, frames)
#
**Arguments**
#
particles	          A ParticleSet object
#
L	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Maximum number of pixels an object can move in two consecutive frames
#
R	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Linkrange, i.e. the number of consecutive frames to search for potential candidate links
#
epsilon1	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;A numeric value, to be used in the formula. Jitter for allowing angular displacements
#
epsilon2	 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp A numeric value, to be used in the formula. Jitter for allowing spatial displacements
#
lambda1	   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A numeric value. Multiplicative factor for the penalty function
#
lambda2	   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A numeric value. Multiplicative factor applied to the angular displacement
#
penaltyFunction &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; A function structured in such a way to be applied as penalty function in the linking
#
include.area	 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     Logical, whether to include also area change of the particles in the cost function calculation
#
frames	     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;       The Frames object that the ParticleSet object is derived from
#
verboseOutput	  &nbsp;&nbsp;&nbsp;    Logical, whether the output should report additional intermediate steps. For debugging use mainly.
#
prog	     &nbsp;&nbsp;&nbsp;&nbsp;         Logical, whether the a progress bar should be shown during the tracking phase
#
include.intensity	 &nbsp;&nbsp;&nbsp;&nbsp; Logical, whether to include also intensity change of the particles in the cost function calculation
#
**Value**
#
A TrajectorySet object
#
**Author**
#
Ayan Bandyopadhyay, Bellarmine College Prep
#
**Examples**
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
**Description**
#
Green pixels, or pixels with H value between 61/360 and 140/360 and S value over 0.15 in the HSV color scheme, are treated as foreground. All other pixels are treated as background. This is useful in detecting the Green Flourescence Protein used to denote cell viability.
#
**Usage**
#
greenThresh(frames)
#
**Arguments**
#
frames	             A Frames object with all 3 color channels(R,G,B)
#
**Value**
#
**Author**
#
**Examples**


###trajMeasures
**Description**
**Usage**

**Arguments**
**Value**
**Author**
**Examples**
