# pso - a simple particle swarm optimization 

This is a basic particle swarm optimization I wrote as an exercise for university.
It is implemented in simple and clean standard C.
Included are 4 different optimization problems: sphere, rosenbrock, rastrigin, schwefel.

Particles leaving the search space are evaluated to positive infinity.

## Installation

Clone the repository.
`git clone https://muhq.space/software/pso.git`

Build pso by running `make` in the repository.

## Usage

	Usage: ./pso [OPTIONS]
	OPTIONS:
		 -p number of particles
		 -h print this help and exit
		 -t number of threads
		 -i iterations
		 -f the function to optimize
		     0 - sphere function
		     1 - rosenbrock function
		     2 - rastrigin function
		     3 - schwefel function
		 -a the a parameter
		 -b the b parameter
