# To make fish only pursue profitable maneuvers

* Gill raker limitations -- just exclude prey types beyond gill raker limits straight-up.

* Note: we have fish just not detect prey if it would be unprofitable to pursue on a per-encounter basis. That's the easiest place to build it into the model in a single spot without messing up anything else.

* Is there any way to make visual acuity reflect variation in the size of items in a prey class? So some proportion of them is detectable, based on whatever proportion of them would have a given angular size, without having abrupt transitions?

# Remaining additions

  
* Maybe try binary masks for maneuvers and visual acuity that I can use to limit foraging and also turn them on and off to see how limiting those factors are?
  
# Other 

* Make sure velocity can be specified by user rather than optimized. Ideally, make this generalizable
to optimize and not optimize any combination of variables.
* Revisit customizing cache spatial resolution to fish size

* Maybe make maneuver cost calls interpolate between the two nearest velocities rather than rounding to the nearest one -- although a max difference of 2.5 cm/s really isn't all that important.

# Things to check

# Printouts

* Make a printout of an optimized maneuver say if it's hitting bounds.

## Minor speed tweaks maybe

* Have the integrands check for detection probability 0, energy value 0, etc, and skip the expensive evaluations in those cases

## Convenience / low priority / when distracted

* Make getters and arrange public/private
* Move DIAG variables into compile-time macros when code is finalized

# Links

* C optimization tips https://people.cs.clemson.edu/~dhouse/courses/405/papers/optimize.pdf

# Optimization possibilities

* GAlib, C++ library for standard genetic algorithms, from the Dos/Win3 era http://lancet.mit.edu/galib-2.4/
* Sago, C++ library for differential evolution with adaptive control parameters, from Windows XP era http://www1.osu.cz/home/habibal/optimization/cpp_library.html

* Milsto's Differential Evolution algorithm, C++, single-header library from 2017, but the algorithm isn't that new: https://github.com/milsto/differential-evolution

* Another DE library, full featured and fancier but still maybe without inequality constraints? Used to optimize trading algorithms: http://www.amichel.com/de/doc/html/dd/d53/overview.html

* PSO: particle swarm optimization in C (no inequality constraints apparently) https://github.com/kkentzo/pso

* eEAag, 2010 contest winner for constrained evolutionary algorithms, http://www.ints.info.hiroshima-cu.ac.jp/~takahama/eng/papers/eDEag-CEC2010.pdf, http://www.ints.info.hiroshima-cu.ac.jp/~takahama/eng/index.html (paper shows way too many function evaluations for good results!!)

* various ant colony optimizers in C here: http://iridia.ulb.ac.be/~mdorigo/ACO/aco-code/public-software.html

* homepage for grey wolf optimization http://www.alimirjalili.com/GWO.html

* comparative study of nature-based algorithms from 2017... journal / scihub website down though https://www.researchgate.net/publication/320461544_A_Comparative_Study_on_Recently-Introduced_Nature-Based_Global_Optimization_Methods_in_Complex_Mechanical_System_Design

* evoloPy for if I want to do the optimization in Python https://github.com/7ossam81/EvoloPy

# Final calibration of model parameters

* Multi-objective optimization for pareto-optimal solutions? https://en.wikipedia.org/wiki/Pareto_efficiency

