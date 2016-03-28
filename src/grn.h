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

#ifndef GRN_H
#define GRN_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef GRN_ACTIV
#define GRN_ACTIV tanh
#endif

struct grn_network
{
	float *genome;
	float *matrix;
	float t1; // magnitude of the interaction
	float t2; // decay rate
	int size;
};

unsigned int grn_rand(void);
float grn_randf(void);

void grn_create(struct grn_network *net, int size);
void grn_destroy(struct grn_network *net);
void grn_copy(struct grn_network *dest, const struct grn_network *src);

float grn_matrix_l1(const struct grn_network *net);

void grn_matrix_mul(const struct grn_network *net, float *output, const float *input);
void grn_run(const struct grn_network *net, float *output, int step_count);

void grn_mutate_single(struct grn_network *net, float u1, float u2);
void grn_mutate_proba(struct grn_network *net, float pr1, float pr2);

#ifdef __cplusplus
}
#endif

#endif