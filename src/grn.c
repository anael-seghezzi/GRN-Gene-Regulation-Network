/*======================================================================
 GRN gene regulation network
 version 1.0
 Implemented from the paper:
 The evolution of phenotypic correlations and ‘developmental memory'
------------------------------------------------------------------------
 Copyright (c) 2016 Anael Seghezzi

 This software is provided 'as-is', without any express or implied
 warranty. In no event will the authors be held liable for any damages
 arising from the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software
    in a product, an acknowledgment in the product documentation would
    be appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not
    be misrepresented as being the original software.

 3. This notice may not be removed or altered from any source
    distribution.

========================================================================*/

#include <stdlib.h>
#include <memory.h>
#include <math.h>

#include "grn.h"

#define _GRN_CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

static unsigned int _grn_rz = 362436069;
static unsigned int _grn_rw = 521288629;

unsigned int grn_rand(void)
{
   _grn_rz = 36969 * (_grn_rz & 65535) + (_grn_rz >> 16);
   _grn_rw = 18000 * (_grn_rw & 65535) + (_grn_rw >> 16);
   return (_grn_rz << 16) + _grn_rw;
}

float grn_randf(void)
{
   unsigned int u = grn_rand();
   return (u + 1.0) * 2.328306435454494e-10;
}

void grn_create(struct grn_network *net, int size)
{
	net->genome = (float *)calloc(size, sizeof(float));
	net->matrix = (float *)calloc(size * size, sizeof(float));
	net->t1 = 1.0;
	net->t2 = 0.2;
	net->size = size;
}

void grn_destroy(struct grn_network *net)
{
	free(net->genome);
	free(net->matrix);
}

void grn_copy(struct grn_network *dest, const struct grn_network *src)
{
	memcpy(dest->genome, src->genome, sizeof(float) * src->size);
	memcpy(dest->matrix, src->matrix, sizeof(float) * (src->size * src->size));
	dest->t1 = src->t1;
	dest->t2 = src->t2;
	dest->size = src->size;
}

float grn_matrix_l1(const struct grn_network *net)
{
	float *m = net->matrix;
	float l1 = 0;
	int i, s2 = net->size * net->size;

	for (i = 0; i < s2; i++) {
		l1 += (*m) < 0 ? -(*m) : (*m);
		m++;
	}

	return l1 / s2;
}

void grn_matrix_mul(const struct grn_network *net, float *output, const float *input)
{
	float *m = net->matrix;
	float idecay = 1.0 - net->t2;
	int i;

	for (i = 0; i < net->size; i++) {

		float v = 0;
		int j;

		for (j = 0; j < net->size; j++) {
			v += input[j] * (*m);
			m++;
		}

		output[i] = input[i] * idecay + net->t1 * GRN_ACTIV(v);
	}
}

void grn_run(const struct grn_network *net, float *output, int step_count)
{
	float *buffer;
	float *in, *out;
	int i;
	
	// first call
	grn_matrix_mul(net, output, net->genome);
	
	if (step_count == 1)
		return;

	// reccurent call
	buffer = (float *)malloc(net->size * sizeof(float));

	in = output;
	out = buffer;

	for (i = 2; i < step_count; i++) {

		float *t;

		grn_matrix_mul(net, out, in);

		t = out;
		out = in;
		in = t;
	}

	// last call
	grn_matrix_mul(net, out, in);

	// output
	if (out != output)
		memcpy(output, out, net->size * sizeof(float));

	free(buffer);
}

void grn_mutate_single(struct grn_network *net, float u1, float u2)
{
	int i;
	float v;
	
	i = grn_rand() % net->size;
	v = net->genome[i] + (grn_randf() * 2.0 - 1.0) * u1;
	net->genome[i] = _GRN_CLAMP(v, -1, 1);

	i = grn_rand() % (net->size * net->size);
	v = net->matrix[i] + (grn_randf() * 2.0 - 1.0) * u2;
	net->matrix[i] = _GRN_CLAMP(v, -1, 1);
}

void grn_mutate_proba(struct grn_network *net, float pr1, float pr2)
{
	int i;
	int s2 = net->size * net->size;
	
	for (i = 0; i < net->size; i++) {
		if (grn_randf() < pr1) {
			float v = net->genome[i] + (grn_randf() * 2.0 - 1.0);
			if (v > -1.0 && v < 1.0)
				net->genome[i] = v;
		}
	}

	for (i = 0; i < s2; i++) {
		if (grn_randf() < pr2) {
			float v = net->matrix[i] + (grn_randf() * 2.0 - 1.0);
			if (v > -1.0 && v < 1.0)
				net->matrix[i] = v;
		}
	}
}
