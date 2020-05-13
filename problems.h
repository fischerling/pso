// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright 2020 Florian Fischer
#pragma once

#include "pso.h"

typedef struct problem {
	char* name;
	interval_t interval;
	double(*func)(particle_t*);
} problem_t;

static const interval_t sphere_interval = {-500, 500};
static double sphere_func(particle_t* particle) {
	size_t size = particle->x->size;
	double sum = 0;

	for (size_t i = 0; i < size; ++i) {
		double x = particle->x->elements[i];
		if (x < sphere_interval.start || x > sphere_interval.end)
			return INFINITY;

		sum += pow(x, 2);
	}
	return sum;
}
static const problem_t sphere = {"sphere", sphere_interval, sphere_func};


static const interval_t rosenbrock_interval = {-30, 30};
static double rosenbrock_func(particle_t* particle) {
	size_t size = particle->x->size;
	double sum = 0;

	for (size_t i = 0; i < size - 1; ++i) {
		double x = particle->x->elements[i];
		double next_x = particle->x->elements[i + 1];

		if (x < rosenbrock_interval.start || x > rosenbrock_interval.end
		        || next_x < rosenbrock_interval.start || next_x > rosenbrock_interval.end)
			return INFINITY;

		sum += 100 * pow(next_x - pow(x, 2), 2) + pow((1 - x), 2);
	}
	return sum;
}
static const problem_t rosenbrock = {"rosenbrock", rosenbrock_interval, rosenbrock_func};

static const interval_t rastrigin_interval = {-5.12, 5.12};
static double rastrigin_func(particle_t* particle) {
	size_t size = particle->x->size;
	double sum = 0;

	for (size_t i = 0; i < size; ++i) {
		double x = particle->x->elements[i];

		if (x < rastrigin_interval.start || x > rastrigin_interval.end)
			return INFINITY;

		sum += pow(x, 2) - 10 * cos(2 * M_PI * x);
	}
	return 10 * size + sum;
}
static const problem_t rastrigin = {"rastrigin", rastrigin_interval, rastrigin_func};

static const interval_t schwefel_interval = {-500, 500};
static double schwefel_func(particle_t* particle) {
	size_t size = particle->x->size;
	double sum = 0;

	for (size_t i = 0; i < size; ++i) {
		double x = particle->x->elements[i];

		if (x < schwefel_interval.start || x > schwefel_interval.end)
			return INFINITY;

		sum += -1 * x * sin(sqrt(fabs(x)));
	}
	return sum;
}
static const problem_t schwefel = {"schwefel", schwefel_interval, schwefel_func};

#define NUM_PROBLEMS 4
static problem_t problems[NUM_PROBLEMS] = {sphere, rosenbrock, rastrigin, schwefel};
