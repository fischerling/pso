// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 Florian Fischer
#pragma once

#include <stdlib.h>
#include <string.h>
#include <time.h>

__thread struct drand48_data *random_data;

static void* check_malloc(size_t size) {
	void *mem = malloc(size);
	if (mem == NULL) {
		perror("malloc");
		exit(1);
	}
	return mem;
}

static void init_random() {
	if (random_data == NULL) {
		random_data = check_malloc(sizeof(struct drand48_data));
		srand48_r(time(NULL), random_data);
	}
}

typedef struct interval {
	double start, end;
} interval_t;

typedef struct vector {
	size_t size;
	double elements[];
} vec_t;

static vec_t* new_vec(size_t size) {
	vec_t* vec = malloc(sizeof(size_t) * sizeof(double) * size);
	if (vec == NULL) {
		perror("malloc");
		exit(1);
	}

	vec->size = size;
	return vec;
}

static void print_vec(vec_t* vec) { 
	if (vec->size == 1) { 
		printf("{%g}", vec->elements[0]); 
	} 

	printf("{%g", vec->elements[0]); 
	for (int i = 1; i < vec->size; i++) { 
		printf(", %g", vec->elements[i]); 
	} 

	printf("}"); 
} 

static void copy_vec(vec_t* src, vec_t* dest) {
	assert(src->size == dest->size);
	memcpy(dest->elements, src->elements, sizeof(double) * src->size);
}

static vec_t* rand_vec(size_t size, vec_t* vec) {
	if (vec == NULL) {
		vec = new_vec(size);
	}

	for (size_t i = 0; i < size; ++i) {
		drand48_r(random_data, &vec->elements[i]);
	}
	return vec;
}

static vec_t* rand_vec_in_interval(size_t size, interval_t interval, vec_t* vec) {
	vec = rand_vec(size, vec);

	for (size_t i = 0; i < size; ++i) {
		double r;
		drand48_r(random_data, &r);
		// get uniform sample from standard uniform distribution
		// https://en.wikipedia.org/wiki/Uniform_distribution_(continuous)#Computational_methods
		vec->elements[i] = interval.start + (interval.end - interval.start) * r;
	}
	return vec;
}

static void member_prod(vec_t* v1, vec_t* v2, vec_t* res) {
	assert(v1->size == v2->size && v2->size == res->size);

	for (size_t i = 0; i < v1->size; ++i) {
		res->elements[i] = v1->elements[i] * v2->elements[i];
	}
}

static void member_sum(vec_t* v1, vec_t* v2, vec_t* res) {
	assert(v1->size == v2->size && v2->size == res->size);

	for (size_t i = 0; i < v1->size; ++i) {
		res->elements[i] = v1->elements[i] + v2->elements[i];
	}
}

static void member_diff(vec_t* v1, vec_t* v2, vec_t* res) {
	assert(v1->size == v2->size && v2->size == res->size);

	for (size_t i = 0; i < v1->size; ++i) {
		res->elements[i] = v1->elements[i] - v2->elements[i];
	}
}

static void scalar_prod(vec_t* v1, double scalar, vec_t* res) {
	assert(v1->size == res->size);

	for (size_t i = 0; i < v1->size; ++i) {
		res->elements[i] = v1->elements[i] * scalar;
	}
}

typedef struct particle {
	double p_val;
	vec_t *v, *x, *p, *r_loc, *r_glob;
} particle_t;

static particle_t* new_particle(size_t size, interval_t interval) {
	particle_t* particle = check_malloc(sizeof(particle_t));

	particle->p_val = INFINITY;

	particle->v = new_vec(size);
	memset(&particle->v->elements[0], 0, sizeof(double) * size);

	particle->x = new_vec(size);
	rand_vec_in_interval(size, interval, particle->x);

	particle->p = new_vec(size);
	copy_vec(particle->x, particle->p);

	particle->r_loc = new_vec(size);
	particle->r_glob = new_vec(size);

	return particle;
}

static void destroy_particle(particle_t* particle) {
	free(particle->x);
	free(particle->v);
	free(particle->p);
	free(particle->r_loc);
	free(particle->r_glob);
	free(particle);
}
static void print_particle(particle_t* particle) __attribute__((unused));

static void print_particle(particle_t* particle) {
	printf("Particle at: ");
	print_vec(particle->x);
	printf(" with v: ");
	print_vec(particle->v);
	printf(" and local best %lf at ", particle->p_val);
	print_vec(particle->p);
	printf("\n");
}

static void evaluate_particle(particle_t* particle, double(*evaluation_func)(particle_t*)) {
	double res = evaluation_func(particle);

	if (particle->p_val == INFINITY || res < particle->p_val) {
		// update local optimum
		particle->p_val = res;
		copy_vec(particle->x, particle->p);
	}
}

static void find_min(particle_t** particles, size_t num_particles, vec_t** optimum_location, double* optimum) {
	particle_t* particle = particles[0];
	*optimum = particle->p_val;
	*optimum_location = particle->p;

	for (size_t i = 1; i < num_particles; ++i) {
		particle = particles[i];
		if (particle->p_val < *optimum) {
			*optimum = particle->p_val;
			*optimum_location = particle->p;
		}
	}
}

static void step(double a, double b_loc, double b_glob,
                 particle_t* particle, vec_t* p_glob) {
	assert(particle->v->size == particle->p->size && particle->p->size == p_glob->size);

	vec_t *v = particle->v;
	vec_t *x = particle->x;
	vec_t *p = particle->p;
	vec_t *r_loc = particle->r_loc;
	vec_t *r_glob = particle->r_glob;

	size_t dimensions = x->size;
	// get new random vectors
	rand_vec(dimensions, r_loc);
	rand_vec(dimensions, r_glob);

	// old speed momentum
	scalar_prod(v, a, v);

	// local optimum attraction
	vec_t* local_opt_attraction = new_vec(dimensions);
	// random local attraction part
	scalar_prod(r_loc, b_loc, r_loc);
	member_diff(p, x, local_opt_attraction);
	member_prod(local_opt_attraction, r_loc, local_opt_attraction);

	// global optimum attraction
	vec_t* global_opt_attraction = new_vec(dimensions);
	// random global attraction part
	scalar_prod(r_glob, b_glob, r_glob);
	member_diff(p_glob, x, global_opt_attraction);
	member_prod(global_opt_attraction, r_glob, global_opt_attraction);

	member_sum(v, local_opt_attraction, v);
	free(local_opt_attraction);

	// new speed
	member_sum(v, global_opt_attraction, v);
	free(global_opt_attraction);

	// new location
	member_sum(x, v, x);
}
