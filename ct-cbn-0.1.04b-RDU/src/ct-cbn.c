
/*
  ct-cbn.c -- Continuous time conjunctive Bayesian networks
  Algorithms for model selection and ML parameter estimation

  Version : 0.1.00
  Author  : Niko Beerenwinkel and Moritz Gerstung
  
  See (N. Beerenwinkel and S. Sullivant, 2007) for further details.
  
  
  Copyright (C)  2007-2008  
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
  int GPS = 0;
  int mode = LEARN_PARAM;
  int B = 0;  // bootstrap samples
  int BM = 0; // bootstrap mode

  int c = 0;
  while((c = getopt(argc, argv, "M:B:R:S:e:N:f:r:XvhP")) != EOF )
    {
      switch(c)
        {
	case 'M':
	  if ((atoi(optarg) == 0) || (atoi(optarg) == 1))
	    BM = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'B':
	  if (atoi(optarg) >= 0)
	    B = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'r':
	  if (atof(optarg) > 0)
	    seed = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'S':
	  if (atof(optarg) > 0)
	    S = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'e':
	  e_flag++;
	  if (atof(optarg) < 1.0)
	    eps = atof(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'N':
	  if (atoi(optarg) >= 0)
	    N_draw = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'R':
	  if (atoi(optarg) > 0)
	    R = atoi(optarg);
	  else
	    error_flag++;
	  break;
	  
	case 'f':
	  filestem = optarg;
	  f_flag++;
	  break;
	  
	case 'X':
	  write_expected_times = 1;
	  break;
	  
	case 'P':
	  GPS = 1;
	  break;
	  
	case 'v':
	  verbose = 1;
	  break;
	  
	case 'h':
	  printf("\nCT-CBN :: Continuous Time Conjunctive Bayesian Networks\n");
	  printf("\n-- Version 0.1.00, Sep 2008\n");

	  printf("\nUsage:\n");
	  printf("  ct-cbn [OPTIONS]\n");

	  printf("\nOptions:\n");
	  printf("  -N <# of samples = 0>\n");
	  printf("     If > 0, the number of samples to draw from the model\n");
	  printf("     If zero (default), the model will be learned from data\n");
	  printf("  -e <eps = 0.0>\n");
	  printf("     If in [0,1], the fraction of violations allowed per edge\n");
	  printf("     If negative, the interval [0,0.5] will be sampled\n"); 
	  printf("       equidistantly with N points; the resulting posets are\n");
	  printf("       written to the _path_ <filestem>\n");
	  printf("     If omitted, the poset in filestem.poset will be used\n"); 
	  printf("       rather than learning the poset from the data\n"); 
	  printf("  -B <# of bootstrap samples>\n");
	  printf("      Number of bootstrap samples (requires eps > 0, N = 0)\n");
	  printf("  -M <bootstrap mode = 0>\n");
	  printf("      0 = re-estimate poset (default)\n");
	  printf("      1 = re-estimate parameters for a fixed poset\n");
	  printf("  -R <# of EM runs = 1)\n");
	  printf("  -S <sampling rate \\lambda_s = 1.0>\n");
	  printf("  -f <filestem>\n");
	  printf("     The directory 'filestem' is required!\n");
	  printf("  -r <random number generator seed>\n");
	  printf("  -v\n");
	  printf("     Verbose\n");
	  printf("  -h\n");
	  printf("     This help\n");
	  
	  printf("\nPrerequisites:\n");
	  printf("  filestem/       - Subdirectory 'filestem'\n");
	  printf("  filestem.pat    - Mutational patterns (genotypes), unless N > 0\n");
	  printf("  filestem.poset  - Event poset used if -e is _not_ set;\n");
	  printf("                    if -e is set, the file is used for determining\n");
	  printf("                    the number of events as specified in the first row\n");
	  printf("  filestem.lambda - Model parameters, if N > 0\n");
	  
	  
	  printf("\nExamples:\n");
	  printf("  ct-cbn -h\n");
	  printf("    Print this message\n");
	  printf("  ct-cbn -f foo\n");
	  printf("    Estimate parameters for poset foo.poset from data foo.pat\n");
	  printf("  ct-cbn -f foo -e 0.05\n");
	  printf("    Estimate poset and parameter from foo.pat\n");
	  printf("  ct-cbn -f foo -e -1\n");
	  printf("    Estimate poset and parameter from foo.pat for a range of eps values\n");
	  printf("  ct-cbn -f foo -e 0.05 -B 1000 -M 0\n");
	  printf("    Draw 1000 bootstrap samples and re-estimate posets\n");
	  printf("  ct-cbn -f foo -B 1000 -M 1\n");
	  printf("    Draw 1000 bootstrap samples and re-estimate parameters of foo.poset\n");
	  printf("  ct-cbn -f foo -N 200\n");
	  printf("    Draw 200 samples from foo.poset and foo.lambda\n");
	  printf("  ct-cbn -f foo -P\n");
	  printf("    Calculate genetic progression score (GPS) for each observation in foo.pat from foo.poset and foo.lambda\n");
	  printf("\n");
	  exit(0);
	  
	default :
	  exit(1);
        }
    }
  

  // bootstrapping does not go with eps<0 nor drawing samples;
  // poset bootstrapping requires "-e <eps>" with eps>0
  if (B > 0)
    if ((eps < 0.0) || (N_draw > 0) || ((BM == LEARN_POSET) && (! e_flag)))
      error_flag++;
  
  if (error_flag)
    {
      printf("Error: Bad parameter values!\n");
      exit(1);
    }
  if (! f_flag)
    {
      printf("Error: No input file specified!\n");
      exit(1);
    }
  
  if (B > 0)
    {
      mode = BM;
    }
  else
    {
      if (e_flag)
	mode = LEARN_BOTH;
      else
	mode = LEARN_PARAM;
    }
  printf("mode = %d,   %d\n", mode, e_flag);
  
  srand(seed);
  RNG = gsl_rng_alloc (gsl_rng_taus);  // global variable
  gsl_rng_set(RNG, seed);  // seed rng

  int i, j, k;
  
  model M;
  read_poset(filestem, &M);

  // precompute binary expansions
  precompute_binary(M.n+1);
  
  M.lin_ext = get_int_array(M.n);  // a linear extension of the poset
  double* lambda = get_double_array(M.n+1);  // Exp rates
  lambda[0] = S;
  
  if (GPS) // Calculate genetic progression score (GPS)
    {
      // Read data
      int N, N_u;
      int** pat = read_patterns(filestem, &N, M.n);
      int* pat_idx = get_int_array(N);
    
      data* D = make_data_set(pat, N, M.n, &N_u, pat_idx);
   
      for (k=0; k<N; k++)
	free(pat[k]);
      free(pat);

      read_lambda(filestem, lambda, M.n);

      // Inititalize data structures
      M.J_P = bfs_order_ideals(M.P, M.n+1, &(M.m), M.lin_ext);
      parents(&M); 
      children(&M);
      construct_sublattices(D, N_u, &M); 
      compatibility(D, N_u, &M);
      
      double* lambda_exit = get_double_array(M.m);
      compute_lambda_exit(lambda, &M, lambda_exit);

      double* Prob = get_double_array(M.m);
      double* Exp = get_double_array(M.m);

      int c=0;
      double gps;
      
      // Calculate GPS for all observations and print to std
      for (k=0; k<N; k++)
	{
	  c = pat_idx[k];
	   if (D[c].is_compatible)
	    {
	      compute_prob(lambda, &M, D[c].J_Q, lambda_exit, Prob); 
	      gps = compute_GPS (lambda, &M, D[c].J_Q, lambda_exit, Prob, Exp);
	      for(i=0;i<M.n;i++)
		printf("%i ",D[c].g[i]);
	      printf(" %f\n", gps);
	    }
	  else
	    {
	      for(i=0;i<M.n;i++)
		printf("%i ",D[c].g[i]);
	      printf(" NA\n");
	    }
	 
	}
    }
  else
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
      
      printf("Poset\tEps\tAlpha\tLoglike\tlambda_s");
      for (i=1; i<=M.n; i++)
	printf("\tlambda_%d", i);
      printf("\n");
      
      if (eps >= 0.0)  // fixed epsilon
	{
	  int b;
	  eps = MIN(eps, 1.0);
	  
	  // single run:
	  select_poset(0, eps, &M, lambda, D, N_u, R, mode,1);
	  if (e_flag) 
	    write_poset(0, filestem, M.P, M.n, -1);
	  write_lambda(filestem, lambda, M.n);
	  
	  // bootstrap runs:
	  double *p_orig = get_double_array(N_u);  // frequencies of original data
	  for (k=0; k<N_u; k++)
	    p_orig[k] = (double) D[k].count;
	  
	  int** bootstrap_count = get_int_matrix(M.n+1, M.n+1);  // all relations
	  int** bootstrap_cover_count = get_int_matrix(M.n+1, M.n+1);  // cover relations only
	  for (i=0; i<=M.n; i++)
	    for (j=0; j<=M.n; j++)
	      bootstrap_count[i][j] = bootstrap_cover_count[i][j] = 0;
	  
	  int** T = get_int_matrix(M.n+1, M.n+1);
	  for (b=1; b<=B; b++)
	    {
	      resample(D, p_orig, N_u);
	      select_poset(b, eps, &M, lambda, D, N_u, R, mode,1);  // e_flag

	      int_matrix_sum(bootstrap_cover_count, M.P, bootstrap_cover_count, M.n+1);
	      transitive_closure(M.P, T, M.n+1);
	      int_matrix_sum(bootstrap_count, T, bootstrap_count, M.n+1);
	      
	      if (e_flag) 
		write_poset(0, filestem, M.P, M.n, b);
	    }	  
	  
	  if ((B > 0) && (BM == 0))
	    {
	      printf("\nbootstrap counts: matrix entry (i,j) counts edge i-->j\n");
	      printf("\nall poset relations =\n");
	      print_int_matrix(bootstrap_count, M.n+1, M.n+1);
	      printf("\ncover relations =\n");
	      print_int_matrix(bootstrap_cover_count, M.n+1, M.n+1);
	    }
	      
	  for (i=0; i<=M.n; i++)
	    {
	      free(T[i]);
	      free(bootstrap_count[i]);
	    }
	  free(T);
	  free(bootstrap_count);
	  free(p_orig);
	}
      else  // sample epsilon
	{
	  double epsilon;
	  int k;
	  for (k=0; k<N; k++)
	    {
	      epsilon = (double) k / (2.0 * (double) N);
	      select_poset(k, epsilon, &M, lambda, D, N_u, R, LEARN_BOTH,1);
	      if (e_flag) 
		write_poset(k, filestem, M.P, M.n, -1);
	    }
	}
      free_data(D, N_u, M.n);
      free_poset(&M);
    }
  else  // sample from given model:
    {
      // construct model:
      M.lin_ext = get_int_array(M.n);  // a linear extension of the poset
      M.J_P = bfs_order_ideals(M.P, M.n+1, &(M.m), M.lin_ext);
      parents(&M);
      children(&M);
      if (verbose)  print_model(&M);
      
      int** pat_draw = get_int_matrix(N_draw, M.n+1);
      double** t_draw = get_double_matrix(N_draw, M.n+1);
      
      read_lambda(filestem, lambda, M.n);
      
      draw_samples(M.P, lambda, M.lin_ext, M.n, pat_draw, t_draw, N_draw);
      write_patterns(filestem, pat_draw, N_draw, M.n);
      write_times(filestem, t_draw, N_draw, M.n);
      
      for (k=0; k<N_draw; k++)
	{
	  free(pat_draw[k]);
	  free(t_draw[k]);
	}
      free(pat_draw);
      free(t_draw);
      
      free_model(&M);
    }

    }
  
  free(lambda);
  
  return 0;    
}



