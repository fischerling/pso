// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 Florian Fischer
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "pso.h"
#include "problems.h"

// defaults
static size_t problem = 0;
static size_t dimensions = 2;
static size_t threads = 1;
static size_t num_particles = 100;
static size_t iterations = 100;
static double a = 0.72984;
static double b_loc = 1.496172;
static double b_glob = 1.496172;

static particle_t** particles;
static vec_t* p_glob;
static double p_glob_val;

static pthread_barrier_t iteration_barrier;
static pthread_barrier_t p_glob_barrier;

static void* thread_func(void* arg) {
	size_t tid = (size_t)arg;

	size_t particles_per_thread = num_particles / threads;
	size_t first_particle = tid * particles_per_thread;
	size_t last_particle = first_particle + particles_per_thread;

	init_random();

	// iteration loop
	for (size_t iteration = 0; iteration < iterations; ++iteration) {
		// particle loop
		for (size_t i = first_particle; i < last_particle; ++i) {
			step(a, b_loc, b_glob, particles[i], p_glob);
			evaluate_particle(particles[i], problems[problem].func);
		}

		pthread_barrier_wait(&iteration_barrier);
		if (tid == 0) {
			// update p_glob
			vec_t* optimum_location;
			find_min(particles, num_particles, &optimum_location, &p_glob_val);
			copy_vec(optimum_location, p_glob);

			// reset barrier
			pthread_barrier_destroy(&iteration_barrier);
			pthread_barrier_init(&iteration_barrier, NULL, threads);
		}

		// wait till p_glob is updated
		pthread_barrier_wait(&p_glob_barrier);
		if (tid == 0) {
			// reset barrier
			pthread_barrier_destroy(&p_glob_barrier);
			pthread_barrier_init(&p_glob_barrier, NULL, threads);
		}
	}

	return NULL;
}

int main(int argc, char* argv[]) {
	// TODO: read parameters from command line
	for(size_t i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-p") == 0) {
			num_particles = atol(argv[i + 1]);
		} else if (strcmp(argv[i], "-i") == 0) {
			iterations = atol(argv[i + 1]);
		} else if (strcmp(argv[i], "-t") == 0) {
			threads = atol(argv[i + 1]);
		} else if (strcmp(argv[i], "-f") == 0) {
			problem = atol(argv[i + 1]);
		} else if (strcmp(argv[i], "-b") == 0) {
			b_loc = strtod(argv[i + 1], NULL);
			b_glob = b_loc;
		} else if (strcmp(argv[i], "-a") == 0) {
			a = strtod(argv[i + 1], NULL);
		} else if (strcmp(argv[i], "-h") == 0) {
			printf("Usage: %s [OPTIONS]\n", argv[0]);
			printf("OPTIONS:\n");
			printf("\t -p number of particles\n");
			printf("\t -h print this help and exit\n");
			printf("\t -t number of threads\n");
			printf("\t -i iterations\n");
			printf("\t -f the function to optimize\n");
			printf("\t     0 - sphere function\n");
			printf("\t     1 - rosenbrock function\n");
			printf("\t     2 - rastrigin function\n");
			printf("\t     3 - schwefel function\n");
			printf("\t -a the a parameter\n");
			printf("\t -b the b parameter\n");
			exit(1);
		}
	}

	assert(num_particles % threads == 0);
	assert(problem < NUM_PROBLEMS);

	printf("threads: %ld, particles: %ld, iterations: %ld, function: %s\n",
	       threads, num_particles, iterations, problems[problem].name);
	printf("a: %lf, b_loc: %lf b_glob: %lf\n",
	       a, b_loc, b_glob);

	pthread_t *threads_array = check_malloc(sizeof(pthread_t*) * threads);

	if (pthread_barrier_init(&iteration_barrier, NULL, threads) != 0) {
		perror("pthread_barrier_init");
		exit(1);
	}

	if (pthread_barrier_init(&p_glob_barrier, NULL, threads) != 0) {
		perror("pthread_barrier_init");
		exit(1);
	}

	init_random();
	particles = check_malloc(sizeof(particle_t*) * num_particles);
	for (size_t i = 0; i < num_particles; ++i) {
		particles[i] = new_particle(dimensions, problems[problem].interval);
	}

	vec_t* optimum_location;
	find_min(particles, num_particles, &optimum_location, &p_glob_val);
	p_glob = new_vec(dimensions);
	copy_vec(optimum_location, p_glob);

	for (size_t i = 0; i < threads; ++i) {
		if (pthread_create(&threads_array[i], NULL, thread_func, (void*)i) != 0) {
			perror("pthread_create");
			exit(1);
		}
	}

	for (size_t i = 0; i < threads; ++i) {
		if (pthread_join(threads_array[i], NULL) != 0) {
			perror("pthread_join");
			exit(1);
		}
	}

	printf("Found optimum %g at ", p_glob_val);
	print_vec(p_glob);
	printf(" after step %ld\n", iterations);

	for (size_t i = 0; i < num_particles; ++i) {
		destroy_particle(particles[i]);
	}

	free(particles);
	free(threads_array);
	return 0;
}
