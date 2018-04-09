# Bugs

* Am I mixing angular area (described in manuscript) with angular width, with respect to either angular resolution or the angular size equation?

# Functionality to add

* Make sure velocity and other variables can be fixed by the user rather than optimized. Right now I just optimize all-or-nothing.

# Possibilities to consider

* Is there any way to make visual acuity reflect variation in the size of items in a prey class? So some proportion of them is detectable, based on whatever proportion of them would have a given angular size, without having abrupt transitions?

## Convenience / low priority / when distracted

* Make more getters and arrange public/private classifications more carefully

# Notes on things to add to manuscript

* We have fish just not detect prey if it would be unprofitable to pursue on a per-encounter basis. That's the easiest place to build it into the model in a single spot without messing up anything else.

* Added gill raker / mouth gape constraints as per previous drift feeding models, except prey types aren't truncated by constraints but fully included/excluded based on their average size. This assumes prey types are divided finely enough that the lack of truncation doesn't represent a major rounding error.

* Mention tau_0 as an index of the fish's aptitude and crypticity as a multiplier on it.

* Write up the angular resolution constraint.

# Useful links

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

