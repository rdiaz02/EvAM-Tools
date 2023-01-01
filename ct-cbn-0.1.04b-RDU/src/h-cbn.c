
/*
  h-cbn.c -- Hidden conjunctive Bayesian networks
  Algorithms for model selection and ML parameter estimation

  Version : 0.1.04
  Author  : Niko Beerenwinkel and Moritz Gerstung

  See  N. Beerenwinkel and S. Sullivant, Biometrika 96(3), 2009, 
          DOI: 10.1093/biomet/asp023, 
  and  M. Gerstung et al., Bioinformatics, 2009
          DOI: 10.1093/bioinformatics/btp505
  for further details.


  Copyright (C)  2007-2011  
  Niko Beerenwinkel                 Moritz Gerstung
  niko.beerenwinkel@bsse.ethz.ch    moritz.gerstung@bsse.ethz.ch

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include "ct-cbn.h"

int main(int argc, char **argv)
{

	// defaults:
	double eps = 0.0;  // e
	int R = 1;  // # of EM runs
	double S = 1.0;  // sampling rate \lambda_s
	char *filestem = "";
	unsigned int seed = (unsigned) time(NULL);  // r, random seed
	verbose = 0;
	int N_draw = 0;  // # of samples to draw
	int error_flag = 0;
	int f_flag = 0;
	int e_flag = 0;
	int gps_flag = 0;
	int l_flag = 0;
	int s_flag = 0;
	int t_flag = 0;
	int n_flag = 0;
	int m_flag = 0;
	int w_flag = 0;
	double T = 0.0;
	int N_iter = 0;
	int only_falsepos = 0;
	int p_flag = 0;
	int mode = LEARN_PARAM;
	double epsilon = 0.0;
	int c = 0;

	while((c = getopt(argc, argv, "e:f:gvhlsFT:N:wpm")) != EOF )
	{
		switch(c)
		{

		case 'e':
			e_flag++;
			if (atof(optarg) < 1.0)
				eps = atof(optarg);
			else
				error_flag++;
			break;

		case 'f':
			filestem = optarg;
			f_flag++;
			break;

		case 'g':
			gps_flag++;
			break;

		case 'l':
			l_flag = 1;
			break;

		case 's':
			s_flag = 1;

		case 'F':
			only_falsepos = 1;
			break;

		case 'T':
			t_flag = 1;
			if (atof(optarg) > 0.0)
				T = atof(optarg);
			else
				error_flag++;
			break;

		case 'N':
			n_flag = 1;
			if (atof(optarg) > 0.0)
				N_iter = atof(optarg);
			else
				error_flag++;
			break;

		case 'v':
			verbose = 1;
			break;

		case 'w':
			w_flag = 1;
			break;

		case 'p':
			p_flag = 1;
			break;

		case 'm':
			m_flag = 1;
			break;

		case 'h':
			printf("\nH-CBN :: Hidden Continuous Time Conjunctive Bayesian Networks\n");
			printf("\n-- Version 0.1.04, Sep 2011\n");

			printf("\nUsage:\n");
			printf("  h-cbn [OPTIONS]\n");

			printf("\nOptions:\n");
			printf("  -f <filestem>\n");
			printf("     The directory 'filestem' is required!\n");
			printf("  -s \n");
			printf("     If set, performes a simulated annealing run starting from the poset <filestem>.poset\n");
			printf("  -T <temp>\n");
			printf("     Sets temperature T of simulated annealing. Default: T=1.\n");
			printf("  -N <steps>\n");
			printf("     Sets the number of simulated annealing steps. Default: N=(N_mutations)^2.\n");
			printf("  -h\n");
			printf("     This help\n");
			printf("\nAdvanced Options:\n");
			printf("  -e <eps = 0.0>\n");
			printf("     Specify value of eps for ct-cbn model selection. Can be useful to determine a starting poset for the h-cbn structure search.\n");
			printf("  -g \n");
			printf("     Compute genetic progression score GPS and MAP estimates. These are written to <filestem.gps> and <filestem.map>, respectively.\n");
			printf("  -l \n");
			printf("     If set, the local log-likelihood is read from the file <filestem>.log. Useful for continuing a simulated annealing run.\n");
			printf("  -m\n     Print the most likely path.\n");
			printf("  -v\n");
			printf("     Verbose\n");
			printf("  -w\n");
			printf("     If set, write estimates of lambda to <filestem>.lambda, the estimated poset to <filestem/00000.poset>, and, if the structure search is activated, the local log-likelihood to <filestem>.log. Repeated after each iteration of the simulated annealing.\n");

			printf("\nPrerequisites:\n");
			printf("  filestem/       - Subdirectory 'filestem'\n");
			printf("  filestem.pat    - Mutational patterns (genotypes)\n");
			printf("  filestem.poset  - Event poset used if -e is _not_ set;\n");
			printf("                    if -e is set, the file is used for determining\n");
			printf("                    the number of events as specified in the first row\n");


			printf("\nExamples:\n");
			printf("  h-cbn -h\n");
			printf("    Print this message\n");
			printf("  h-cbn -f foo\n");
			printf("    Estimate parameters for poset foo.poset from data foo.pat\n");
			printf("  h-cbn -f foo -e 0.05\n");
			printf("    Estimate poset according to ct-cbn and estimated parameters with h-cbn from foo.pat\n");
			printf("  h-cbn -f foo -s -T 10 -N 200");
			printf("    Do a simulated annealing search starting from the poset foo.poset with initial temperature T=10 and do 200 steps.\n");
			printf("  h-cbn -f foo -g ");
			printf("    Estimate parameters and compute the GPS and MAP estimates for foo.poset.\n");

			printf("\n");
			printf("Note:\n  The number of threads can be set via the environment variable OMP_NUM_THREADS.\n\n");

			exit(0);

		default :
			exit(1);
		}
	}

	if(!f_flag) {
	  printf("Error! -f flag must be specified. See h-cbn -h.\n");
	  exit(1);
	}
	srand(seed);
	RNG = gsl_rng_alloc (gsl_rng_taus);  // global variable
	gsl_rng_set(RNG, seed);  // seed rng

	int i, k;

	model M;
	read_poset(filestem, &M);

	// precompute binary expansions
	precompute_binary(M.n+1);

	M.lin_ext = get_int_array(M.n);  // a linear extension of the poset
	double* lambda = get_double_array(M.n+1);  // Exp rates
	lambda[0] = S;

	double total_loglik = 0.0;

	{

		if (N_draw == 0)  // learn model
		{
			int N, N_u;
			int** pat = read_patterns(filestem, &N, M.n);
			int* pat_idx = get_int_array(N);
			data* D = make_data_set(pat, N, M.n, &N_u, pat_idx);
			for (k=0; k<N; k++)
				free(pat[k]);
			free(pat);

			printf("Poset\tEps\tAlpha\tLoglik\tlambda_s");
			for (i=1; i<=M.n; i++)
				printf("\tlambda_%d", i);
			printf("\n");

			if (eps >= 0.0)  // fixed epsilon
			{
				eps = MIN(eps, 1.0);

				/*  single run: */
				select_poset(0, eps, &M, lambda, D, N_u, R, mode,0);

				if(eps > 0){
					read_lambda(filestem, lambda, M.n);
				}

				/* Compute variables */
				double* Prob = get_double_array(M.m);
				double** condprob = get_double_matrix(M.m, N_u);
				double* lambda_exit = get_double_array(M.m);
				compute_lambda_exit(lambda, &M, lambda_exit);
				int* lattice_index = get_int_array(pow2(M.n+1));
				for (i=0; i < M.m; i++)
					lattice_index[M.J_P[i]]=i;
				compute_all_prob(lambda, &M, lambda_exit, Prob, lattice_index);

				/* Compute prob of all observations for ct-cbn;
	   needed to be done here, befor params are re-estimated. */
				if (p_flag){
					double* ProbY_ctcbn = compute_ProbY_ctcbn(&M, Prob);
					write_double_array(filestem, ".prY_ctcbn", ProbY_ctcbn, pow2(M.n));
				}


				if (eps == 0){ // ie, epsilon not specified
					epsilon=0.5; // initial value

					/* Estimate epsilon | lambda */
					if(verbose) printf("\n===EM_epsilon===\n");
					total_loglik = EM_epsilon(&M, D, N_u, Prob, condprob, &epsilon);
					if(verbose) {
						printf("\t%f\t%g\t%g\t", epsilon, pow(1-epsilon,M.n), total_loglik);
						print_double_array(lambda, M.n+1);
					}

					/* Estimate epsilon and lambda */
					int N_compatible = 0;
					int N = 0;
					double alpha;
					for (k=0; k<N_u; k++)
					{
						N += D[k].count;
						N_compatible += D[k].is_compatible * D[k].count;
					}
					alpha = (double) N_compatible / (double) N;

					if(verbose) printf("\n+++EM_EM+++\n");
					epsilon = 0.00001;
					total_loglik = EM_EM(&M, D, N_u, lambda, &epsilon);
					printf("\t%f\t%g\t%g\t",  epsilon, alpha, total_loglik);
					print_double_array(lambda, M.n+1);
				}
				else{ //ie, given epsilon
					epsilon = eps;
					printf("\n+++HM-CBN+++\n");
					total_loglik = compute_total_loglik(&M, D, N_u, Prob, condprob, &epsilon);
					printf("\t%f\t%g\t%g\t", epsilon, pow(1-epsilon,M.n), total_loglik);
					print_double_array(lambda, M.n+1);
				}

				/*Write output?*/
				if (w_flag){
					write_poset(0, filestem, M.P, M.n, -1);
					write_lambda(filestem, lambda, M.n);
				}

				/* Output probability of each observation Y */
				if (p_flag){
					double* ProbY = compute_ProbY(&M, Prob, only_falsepos, epsilon);
					write_double_array(filestem, ".prY", ProbY, pow2(M.n));
				}


				/* GPS */
				if(gps_flag)
				{
				  //printf("\n+++GPS+++\n");

					double* all_GPS = get_double_array(M.m);
					double* cond_GPS = get_double_array(N_u);
					double* loglik = get_double_array(N_u);

					GPS(&M, D, N_u, lambda, epsilon, all_GPS, cond_GPS);
					compute_loglik(&M, D, N_u, Prob, condprob, &epsilon, loglik);

					int* ml_pat = get_int_array(N_u);
					double** exp_pat = get_double_matrix(M.n,N_u);
					compute_hidden_patterns(&M, D, N_u, lambda, epsilon, condprob, ml_pat, exp_pat);

					// Write to file
					write_gps(filestem, cond_GPS, N, pat_idx);

					// Write MAP to file
					char suffix[15] = ".map";
					char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
					strcat(filename, filestem);
					strcat(filename, suffix);

					FILE *output;
					if ( (output = fopen(filename, "w")) == NULL)
					  {
					    fprintf(stderr, "Error:  Could not write to file %s\n", filename);
					    exit(1);
					  }

					for (k=0; k<N;k++)
					{
						c = pat_idx[k]; //Unique observation index
						int ml_gen = M.J_P[ml_pat[c]]; //Lattice entry
						for(i=0;i < M.n + 1;i++)
						  fprintf(output, "%i ", GENOTYPE[ml_gen][i]);
						//printf(" %i ", hamdist(index_of(D[c].g, M.n+1), ml_gen));
						//printf(" ");
						//for(i=0; i < M.n; i++)
						//	printf("%.5f ", exp_pat[i][c]);
						//printf(" %.5f ", cond_GPS[c]);
						//printf(" %.5f", all_GPS[ml_pat[c]]);
						//printf(" %.5f\n", loglik[c]);
						fprintf(output, "\n");

					}
					fclose(output);
				}

				/* Local search */
				if (s_flag)
				{
					if (!n_flag)
						N_iter =  2 * M.n * M.n;
					printf("\n---Local search---\n");
					local_search(&M, D, N_u, lambda, &epsilon, total_loglik, T, N_iter, filestem, l_flag, w_flag);
				}

			}
			free_data(D, N_u, M.n);
			/* Print final poset */
			// for (i = 0; i< pow2(M.n); i++)
			// print_int_array(GENOTYPE[i + pow2(M.n)], M.n+1);
			if(m_flag) ML_path(&M, lambda);
			free_poset(&M);
		}
	}
	free(lambda);
	return 0;
}



