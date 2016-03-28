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

#include <stdio.h>
#include "../src/grn.c"


//#define REGL1 0.1 // L1 regularization


float exp2_fitness(int pattern, float phenotype[8])
{
	float S1[8] = {1, 1, -1, -1, -1, 1, -1, 1};
	float S2[8] = {1, -1, 1, -1, 1, -1, -1, -1};
	float *S;
	float dot = 0;
	int i;

	if (pattern == 0)
		S = S1;
	else
		S = S2;

	for (i = 0; i < 8; i++)
		dot += phenotype[i] * S[i];

	return 1 + dot;
}

int main(int argc, char **argv)
{
	float phenotype[8];
	struct grn_network net, netcpy;
	int pattern = 0;
	int i;


	// create networks
	grn_create(&net, 8);
	grn_create(&netcpy, 8);


	// run generations
	for (i = 0; i < 800000; i++) {

		float fitness, f1, f2;

		if (i > 0 && (i % 2000) == 0) // alternate pattern
			pattern = !pattern;

		grn_copy(&netcpy, &net);
		grn_mutate_single(&net, 0.1, 0.0067);
		grn_mutate_single(&netcpy, 0.1, 0.0067);

		grn_run(&net, phenotype, 10);
		f1 = exp2_fitness(pattern, phenotype);
		
		grn_run(&netcpy, phenotype, 10);
		f2 = exp2_fitness(pattern, phenotype);

#ifdef REGL1
		f1 += (1.0 - grn_matrix_l1(&net)) * REGL1;
		f2 += (1.0 - grn_matrix_l1(&netcpy)) * REGL1;
#endif

		if (f2 > f1) {
			grn_copy(&net, &netcpy);
			fitness = f2;
		}
		else {
			fitness = f1;
		}

		printf("%.5d: fitness = %f\n", i, fitness);
	}


	grn_destroy(&net);
	grn_destroy(&netcpy);
	return EXIT_SUCCESS;
}
