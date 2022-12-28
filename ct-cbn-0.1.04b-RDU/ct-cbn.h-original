#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include <omp.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#include "queue.h"


#define DOUBLE_FORMAT "%.10g"
#define EM_ACCURACY 1e-8
#define EM_MAX_ITER 100000

#define MAX(a,b) ( (a) > (b) ? (a) : (b) ) 
#define MIN(a,b) ( (a) < (b) ? (a) : (b) ) 


enum L_MODE {LEARN_POSET, LEARN_PARAM, LEARN_BOTH};

int verbose;
int write_expected_times;

gsl_rng *RNG;  // random number generator

int** GENOTYPE;  // GENOTYPE[i] is the integer i in binary (as int array)


typedef struct {
	int** P;  // cover relations of the event poset
	int n;  // number of events
	int* lin_ext;  // linear extension of the poset

	int* J_P;  // lattice of order ideals (genotypes)
	int m;  // lattice size

	double eps;  // error tolerance used for constructing P

	int* N_pa;  // number of parents for each genotype
	int** pa;  // list of parents
	int** pa_diff;  // index in which parent and child differ

	int* N_ch; // number of children in lattice
	int** ch; // list of children
	int** ch_diff;  // index in which parent and child differ
} model;


typedef struct {
	int *g;  // genotype
	double* t;  // time of each event
	int is_compatible;  // compatible with model?
	int** Q;  // induced refinement of poset
	int* J_Q;  // corresponding sublattice
	int count;  // number of observations of this type
} data;

typedef struct {
	model* M; //model
	data* D; //data
	int* N_u; //No obs
} gsl_parameters;


int* get_int_array(const int n)
{
	int* x = calloc(n, sizeof(int));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	return x;
}


unsigned int* get_uint_array(const int n)
{
	unsigned int* x = calloc(n, sizeof(unsigned int));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	return x;
}


double* get_double_array(const int n)
{
	double* x = calloc(n, sizeof(double));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	return x;
}


int** get_int_matrix(const int m, const int n)
{
	int** x = malloc(m * sizeof(int *));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	int i;
	for (i=0; i<m; i++)
	{
		x[i] = calloc(n, sizeof(int));
		if (x[i] == NULL)
		{
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
	}

	return x;
}


double** get_double_matrix(const int m, const int n)
{
	double** x = malloc(m * sizeof(double *));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	int i;
	for (i=0; i<m; i++)
	{
		x[i] = calloc(n, sizeof(double));
		if (x[i]  == NULL)
		{
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
	}

	return x;
}

double*** get_double_cube(const int m, const int n, const int l)
{
	double*** x = malloc(m * sizeof(double **));

	if (x == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	int i;
	for (i=0; i<m; i++)
	{
		x[i] = get_double_matrix(n,l);
		if (x[i]  == NULL)
		{
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
	}

	return x;
}

void print_int_array(int* x, int n)
{
	int j;

	if (n > 0)
	{
		for (j=0; j<n-1; j++)
			printf("%d ", x[j]);
		printf("%d", x[n-1]);
	}
	printf("\n");
}



void print_int_matrix(int** X, int m, int n)
{
	int i;

	for (i=0; i<m; i++)
		print_int_array(X[i], n);

}


void print_double_array(double* x, int n)
{
	int j;

	if (n > 0)
	{
		for (j=0; j<n; j++)
		{
			printf(DOUBLE_FORMAT, x[j]);
			if (j < n-1)
				printf("\t");
		}
	}
	printf("\n");
}



void print_double_matrix(double** X, int m, int n)
{
	int i;

	for (i=0; i<m; i++)
		print_double_array(X[i], n);

}



void write_poset(int k, char* filestem, int** P, int n, int b)
{
	int i, j;

	char filename[255];

	if (b >= 0)
		sprintf(filename, "%s/b%05d.poset", filestem, b);
	else
		sprintf(filename, "%s/%05d.poset", filestem, k);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL )
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		fprintf(stderr, "        Make sure the directory '%s' exists.\n", filestem);
		exit(1);
	}

	fprintf(output, "%d\n", n);
	for (i=1; i<=n; i++)
		for (j=1; j<=n; j++)
			if (P[i][j])
				fprintf(output, "%d %d\n", i, j);

	fprintf(output, "0\n");
	fclose(output);
}



void write_patterns(char* filestem, int** pat, int N, int n)
{
	int i, j;

	char suffix[15] = ".sim.pat";
	//char suffix[15] = ".pat";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL )
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	fprintf(output, "%d %d\n", N, n+1);
	for (i=0; i<N; i++)
	{
		for (j=0; j<n; j++)
			fprintf(output, "%d ", pat[i][j]);
		fprintf(output, "%d\n", pat[i][n]);
	}

	fclose(output);
}



void write_times(char* filestem, double** t, int N, int n)
{
	int i, j;

	char suffix[15] = ".sim.time";
	//char suffix[15] = ".time";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL )
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	fprintf(output, "%d %d\n", N, n+1);
	for (i=0; i<N; i++)
	{
		for (j=0; j<n; j++)
			fprintf(output, "%lf ", t[i][j]);
		fprintf(output, "%lf\n", t[i][n]);
	}

	fclose(output);
}





inline void genotype_of(int index, int* x, int n)
{
	int i;

	for (i=n-1; i>=0; i--)
	{
		x[i] = index % 2;
		index = index / 2;
	}

}



int pow2(int k)
{
	int i, pow = 1;

	for (i=0; i<k; i++)
		pow *= 2;

	return(pow);
}



void precompute_binary(const int n) 
{
	int i, j;
	int m = pow2(n);
	int *g = get_int_array(n);

	GENOTYPE = get_int_matrix(m, n);

	for (i=0; i<m; i++)
	{
		genotype_of(i, g, n);
		for (j=0; j<n; j++)
			GENOTYPE[i][j] = g[j];
	}
	free(g);
}



int index_of(int* x, int n)
{
	int i, index = 0;

	for (i=0; i<n; i++)
	{
		if (x[i] == 1)
			index += pow2(n-1-i);
	}

	return index;
}



void free_poset(model* M)
{
	int j;

	for (j=0; j<=M->n; j++)
		free(M->P[j]);
	free(M->P);
	free(M->lin_ext);
}



void free_lattice(model* M)
{
	int i;

	free(M->J_P);
	for (i=0; i<M->m; i++)
	{
		free(M->pa[i]);
		free(M->pa_diff[i]);
		free(M->ch_diff[i]);
	}
	free(M->pa);
	free(M->pa_diff);
	free(M->N_pa);
	free(M->ch_diff);
}



void free_model(model* M)
{
	free_poset(M);
	free_lattice(M);
}



void print_model(model* M)
{
	int i;
	int m = M->m;
	int n = M->n;
	int* g = get_int_array(n+1);

	printf("\n\nMODEL\n\n");
	printf("\nP =\n");
	print_int_matrix(M->P, n+1, n+1);
	printf("\n");

	printf("lattice size, m = %d\nsorted lattice = \n", m);
	for (i=0; i<m; i++)
	{
		printf("%d\t", i);
		genotype_of(M->J_P[i], g, n+1);
		print_int_array(g, n+1);
	}
	//print_int_array(M->J_P, m);

	printf("\nlinear extension of the poset = ");
	print_int_array(M->lin_ext, n);

	for (i=0; i<m; i++)
	{
		printf("\nparents of %5d = ", i);
		//genotype_of(M->J_P[i], g, n+1);
		//print_int_array(g, n+1);
		print_int_array(M->pa[i], M->N_pa[i]);
		printf("differing events = ");
		print_int_array(M->pa_diff[i], M->N_pa[i]);
	}
	printf("\n\n");

	for (i=0; i<m; i++)
	{
		printf("differing events to children of %5d = ", i);
		//genotype_of(M->J_P[i], g, n+1);
		//print_int_array(g, n+1);
		print_int_array(M->ch_diff[i], M->N_ch[i]);
	}

	printf("\n");

	free(g);
}



void print_data(data* D, int N_u, int n, int m)
{
	int k;

	printf("\n\nDATA\n\n");
	for (k=0; k<N_u; k++)
	{
		printf("# %d\n", k);
		printf("g = ");  print_int_array(D[k].g, n+1);
		printf("t = ");  print_double_array(D[k].t, n+1);
		printf("Q =\n");  print_int_matrix(D[k].Q, n+1, n+1);
		printf("J_Q = ");  print_int_array(D[k].J_Q, m);
		printf("count = %d\n", D[k].count);
		printf("is compatible = %d\n", D[k].is_compatible);
		printf("--------------------\n");
	}

}



void free_data(data* D, int N_u, int n)
{
	int k, j;

	for (k=0; k<N_u; k++)
	{
		free(D[k].g);
		free(D[k].t);
		for (j=0; j<=n; j++)
			free(D[k].Q[j]);
		free(D[k].Q);
		free(D[k].J_Q);
	}
	free(D);
}



int** read_patterns(char* filestem, int* N, int n)
{
	int j, k, p;

	char suffix[15] = ".pat";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *input;
	if ( (input = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Error:  Could not read %s\n", filename);
		exit(1);
	}

	/* Read dimensions */
	fscanf(input, "%d %d", N, &p);
	if (verbose) printf("\nreading data from file %s :  %d samples, %d events ...\n\n", filename, *N, p-1);
	if (*N < 1)
	{
		fprintf(stderr, "Error:  Less than one data point!\n");
		exit(1);
	}
	if (n+1 != p)
	{
		fprintf(stderr, "Error:  Number of events in poset and data do not match!\n");
		exit(1);
	}

	int** pat = get_int_matrix(*N, p);

	/* Read patterns */
	int x;
	for (k=0; k<*N; k++)
	{
		for (j=0; j<=n; j++)
		{
			if (fscanf(input,"%d ", &x) == 1)
			{
				if ((x != 0) && (x != 1))
					x = -1;
				pat[k][j] = x;
			}
			else
			{
				fprintf(stderr, "Error reading data from %s!\n", filename);
				exit(1);
			}

		}
	}

	return pat;
}



double** read_times(char* filestem, int* N, int n)
{
	int i, j, p;

	char suffix[15] = ".time";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *input;
	if ( (input = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Error:  Could not read %s\n", filename);
		exit(1);
	}

	/* Read dimensions */
	fscanf(input, "%d %d", N, &p);
	printf("read times: %d samples, %d events\n\n", *N, p-1);
	if (*N < 1)
	{
		fprintf(stderr, "Error:  Less than one data point!\n");
		exit(1);
	}
	if (n+1 != p)
	{
		fprintf(stderr, "Error:  Number of events in poset and data do not match!\n");
		exit(1);
	}

	double** t = get_double_matrix(*N, p);

	/* Read patterns */
	double x;
	for (i=0; i<*N; i++)
	{
		for (j=0; j<=n; j++)
		{
			if (fscanf(input,"%lf ", &x) == 1)
			{
				t[i][j] = x;
			}
			else
			{
				fprintf(stderr, "Error reading data from %s!\n", filename);
				exit(1);
			}

		}
	}

	return t;
}



void read_poset(char* filestem, model* M)
{
	int left, right;

	char suffix[15] = ".poset";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *input;
	if ( (input = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Error:  Could not read %s\n", filename);
		exit(1);
	}

	/* Read number of relations */
	int n;
	fscanf(input, "%d", &n);
	if (verbose)  printf("n = %d events\n\n", n);
	if ((n < 1) || (n > 25))
	{
		fprintf(stderr, "Error:  Number of events is %d.  Supported range is {1, ..., 14}.\n", n);
		exit(1);
	}

	M->n = n;
	M->P = get_int_matrix(n+1, n+1);

	/* Read partial orderings from file */
	fscanf(input,"%d %d", &left, &right);
	while (left != 0)
	{
		if (verbose)  printf("%d --> %d\n", left, right);  // i.e., left < right
		if ((left > n) || (right > n) || (left < 0) || (right < 1))
		{
			fprintf(stderr, "Error:  Undefined event in %s!\n", filename);
			exit(1);
		}
		M->P[left][right] = 1;
		fscanf(input,"%d %d", &left, &right);
	}

	fclose(input);
}



void read_lambda(char* filestem, double* lambda, int n)
{
	int j;

	char suffix[15] = ".lambda";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *input;
	if ( (input = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Error:  Could not read %s\n", filename);
		exit(1);
	}

	double x;
	for (j=0; j<=n; j++)
	{
		fscanf(input, "%lf", &x);
		lambda[j] = x;
	}

	fclose(input);
}



void write_lambda(char* filestem, double* lambda, const int n)
{
	int i;

	char suffix[15] = ".lambda";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	for (i=0; i<=n; i++)
	{
		fprintf(output, "%lf\n", lambda[i]);
	}

	fclose(output);
}

void write_gps(char* filestem, double* cond_GPS, int N, int* pat_idx)
{
	int k,c;

	char suffix[15] = ".gps";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	for (k=0; k<N; k++)
	{
		c = pat_idx[k];
		fprintf(output, "%lf\n", cond_GPS[c]);
	}

	fclose(output);
}

void write_double_array(char* filestem, char suffix[15], double* array, int n)
{
	int i;
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	for (i=0; i<n; i++){
		fprintf(output, DOUBLE_FORMAT, array[i]);
		fprintf(output, "\n");
	}

	fclose(output);

}


void print_genotype(int* x, int n)
{
	int i;

	for (i=0; i<n; i++)
		printf("%d", x[i]);

}



int* bfs_order_ideals(int** poset, const int len, int* count, int* lin_ext)
{
	// implements BFS in the genotype lattice

	int i, j, k, is_compatible;
	int lin_ext_size = 0;

	queue q;
	int g_idx, new_idx;  // current, new genotype index
	int* g = get_int_array(len);  // current genotype

	int* lattice = malloc(sizeof(int));

	int* added = get_int_array(pow2(len));  // records added genotypes
	/* TODO: make this dynamic, e.g., using a list */

	init_queue(&q);
	new_idx = 0;  // wild type 0...0
	enqueue(&q, new_idx);
	added[new_idx] = 1;

	*count = 0;
	while (empty(&q) == FALSE) {
		g_idx = dequeue(&q);
		genotype_of(g_idx, g, len);

		// visit...
		(*count)++;
		if ((lattice = (int*) realloc(lattice, (*count) * sizeof(int))) == NULL)
		{
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
		lattice[(*count)-1] = g_idx;

		// linear extension
		for (i=1; i<len; i++)  // exclude 0
		{
			if (g[i] == 1)
			{
				int is_in = 0;
				for (j=0; j<lin_ext_size; j++)
					if (lin_ext[j] == i)
					{
						is_in = 1;
						break;
					}
				if (! is_in)  // add to linear extension:
				{
					lin_ext[lin_ext_size] = i;
					lin_ext_size++;
				}
			}
		}


		// generate children:
		for (i=0; i<len; i++)
		{
			if (g[i] == 0)
			{
				g[i] = 1;  // try this new event
				new_idx = index_of(g, len);

				// check compatibility:
				is_compatible = 1;
				for (k=0; k<len; k++)
				{
					if ((poset[k][i] == 1) && (g[k] == 0))
					{
						is_compatible = 0;
						break;
					}
				}

				if ((is_compatible) && (! added[new_idx]))
					/* add if compatible and really new */
				{
					enqueue(&q, new_idx);
					added[new_idx] = 1;
				}

				g[i] = 0;  // undo event
			}
		}
	}

	free(g);
	free(added);

	return(lattice);
}



int norm1(const int g_idx, const int n)
{
	int i;
	int *g = get_int_array(n);

	genotype_of(g_idx, g, n);

	int norm1 = 0;
	for(i=0; i<n; i++)
		norm1 += g[i];

	free(g);

	return norm1;
}


double power(double m, int n)
{
	int i;
	double x=1;
	for(i=0;i<n;i++)
		x*=m;
	return x;
}

unsigned hamdist(unsigned x, unsigned y)
{
	unsigned dist = 0, val = x ^ y;

	while(val)
	{
		++dist;
		val &= val - 1;
	}

	return dist;
}

int hamming_distance(int g_idx, int h_idx, int* diff_idx, int n)
{
	int i;
	int *g = get_int_array(n);
	int *h = get_int_array(n);

	genotype_of(g_idx, g, n);
	genotype_of(h_idx, h, n);

	int dist = 0;
	for(i=0; i<n; i++)
		if (g[i] != h[i])
		{
			dist++;
			*diff_idx = i;
		}

	free(g);
	free(h);

	return dist;
}



void parents(model* M)
{
	int i, j, k, c;
	int m = M->m;
	int n = M->n;

	M->N_pa = get_int_array(m);
	M->pa = malloc(m * sizeof(int*));
	M->pa_diff = malloc(m * sizeof(int*));
	if ((M->pa == NULL) || (M->pa_diff == NULL))
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	// count parents:
	for (i=0; i<m; i++)
	{
		M->N_pa[i] = 0;  // number of parents
		k = i-1;
		while ((k >=0) && (norm1(M->J_P[k], n+1) >= (norm1(M->J_P[i], n+1) - 1)))
		{
			//if (hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
			if (hamdist(M->J_P[i], M->J_P[k]) == 1)
			{
				M->N_pa[i]++;  // found a parent!
			}
			k--;
		}
	}

	// list parents:
	for (i=0; i<m; i++)
	{
		M->pa[i] = get_int_array(M->N_pa[i]);
		M->pa_diff[i] = get_int_array(M->N_pa[i]);
		// generate parents in sublattice
		k = i-1;
		c = 0;
		while ((k >=0) && (norm1(M->J_P[k], n+1) >= (norm1(M->J_P[i], n+1) - 1)))
		{
			if(hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
			{
				M->pa[i][c] = k;  // save index of parent in J_P[]
				M->pa_diff[i][c] = j;  // save index by which parent and child differ
				c++;
			}
			k--;
		}
	}

}



void children(model* M)
{
	int i, j, k, c;
	int m = M->m;
	int n = M->n;

	M->N_ch = get_int_array(m);
	M->ch = malloc(m * sizeof(int*));
	M->ch_diff = malloc(m * sizeof(int*));
	if (M->ch_diff == NULL || M->ch == NULL)
	{
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	// count children:
	for (i=0; i<m; i++)
	{
		M->N_ch[i] = 0;
		/* generate children in lattice */
		k = i + 1;
		while ( (k < m) && (norm1(M->J_P[k], n+1) <= (norm1(M->J_P[i], n+1) + 1)) )
		{
			//if (hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
			if (hamdist(M->J_P[i], M->J_P[k]) == 1)
			{
				M->N_ch[i]++;
			}
			k++;
		}
	}

	// list children:
	for (i=0; i<m; i++)
	{
		M->ch[i] = get_int_array(M->N_ch[i]);
		M->ch_diff[i] = get_int_array(M->N_ch[i]);
		/* generate children in lattice */
		k = i + 1;
		c = 0;
		while ( (k < m) && (norm1(M->J_P[k], n+1) <= (norm1(M->J_P[i], n+1) + 1)) )
		{
			if (hamming_distance(M->J_P[i], M->J_P[k], &j, n+1) == 1)
			{
				M->ch[i][c] = k;  // save index of child in J_P[]
				M->ch_diff[i][c] = j;
				c++;
			}
			k++;
		}
	}

}



void compute_prob(double* lambda, model* M, int* sublattice, double* lambda_exit, double* Prob)
{
	int i, j, k, c;

	Prob[0] = 1.0;
	for (i=1; i<M->m; i++)
	{
		Prob[i] = 0.0;
		if (sublattice[i])
		{
			for (c=0; c<M->N_pa[i]; c++)
			{
				k = M->pa[i][c];  // index of c-th parent of i
				if (sublattice[k])
				{
					j = M->pa_diff[i][c];  // index of differing event
					Prob[i] += (lambda[j] / lambda_exit[k]) * Prob[k];
				}
			}
		}
	}

}



void compute_exp(double* lambda, model* M, int* g, int* sublattice, double* lambda_exit, double* Prob, double** Exp)
{
	/*
      Compute  E[X_i - max_{j \in \pa(i)} X_j | X compatible with Q]
      ( in code:  pos = i,  i = index over lattice elements (S) )
	 */

	int i, j, k, c, pred, pos;

	for (pos=0; pos<=M->n; pos++)
	{
		for (i=0; i<M->m; i++)
		{
			Exp[pos][i] = 0.0;
			if (sublattice[i])
			{
				//print_int_array(GENOTYPE[M->J_P[i]],M->n+1);
				for (c=0; c<M->N_pa[i]; c++)
				{
					k = M->pa[i][c];  // index of c-th parent (in J_P) of i
					if (sublattice[k])
					{
						j = M->pa_diff[i][c];  // index of differing event
						g = GENOTYPE[M->J_P[k]];
						int all_pred_in_k = 1; // all predecessor events in J_P[k]?
						pred = 1;
						while (all_pred_in_k && (pred <= M->n))
						{
							all_pred_in_k = ((! M->P[pred][pos]) || (g[pred]));
							pred++;
						}

						Exp[pos][i] += (lambda[j] / lambda_exit[k]) * Exp[pos][k];
						if (all_pred_in_k) {
							Exp[pos][i] += ((1 - g[pos]) * lambda[j] * Prob[k]) /
									(lambda_exit[k] * lambda_exit[k]);
						}
					}
				}
			}
		}
	}
}

// Computes Prob[X], ie the probability of observing the pattern X
void compute_all_prob(double* lambda, model* M, double* lambda_exit, double* Prob, int* lattice_index)
{
	int m = M->m;
	int n = M->n;

	int i,j,k,c;
	int power = pow2(n);

	/*
  int* lattice_index = get_int_array(pow2(n+1));
  for (i=0; i < m; i++)
    lattice_index[M->J_P[i]]=i;
	 */

	// COMPUTE PROB[X,LAMBDA] FOR ALL X IN J(P)

	Prob[0] = 1.0;
	Prob[lattice_index[M->J_P[0] + power]] = Prob[0] * lambda[0] / lambda_exit[0];
	for (i=1; i < m; i++)
	{
		if ( !GENOTYPE[M->J_P[i]][0] )
		{
			Prob[i] = 0.0;
			for (c=0; c < M->N_pa[i]; c++)
			{
				k = M->pa[i][c];  // index of c-th parent of i
				if ( !GENOTYPE[M->J_P[k]][0] )
				{
					j = M->pa_diff[i][c];  // index of differing event
					Prob[i] += (lambda[j] / lambda_exit[k]) * Prob[k];
				}
			}
			// Compute probability to OBSERVE pattern i.
			Prob[lattice_index[M->J_P[i] + power]] = Prob[i] * lambda[0] / lambda_exit[i];
		}
	}
}

//Check if l is >= i
int is_after(model* M, int l, int i)
{
	int a;
	if (M->J_P[l] < M->J_P[i])
		return 0;
	else
		for ( a=0; a<=M->n; a++ )
			if ( GENOTYPE[M->J_P[i]][a] && !GENOTYPE[M->J_P[l]][a] )
				return 0;
	return 1;
}

//Check if l > i
int is_strict_after(model* M, int l, int i)
{
	int a;
	if (M->J_P[l] <= M->J_P[i])
		return 0;
	else
		for ( a=0; a<=M->n; a++ )
			if ( GENOTYPE[M->J_P[i]][a] && !GENOTYPE[M->J_P[l]][a] )
				return 0;
	return 1;
}


// Compute Prob[G|X], the probability of all Genotypes given a censored observation X
void compute_censored_prob(double* lambda, model* M, double* lambda_exit, double* Prob, double** censprob, int* lattice_index)
{
	int m = M->m;
	//int n = M->n;

	int i,j,k,c,l;

	for (i=0; i < m; i++)
	{
		for (l = 0; l < m; l++) // Clear everything
			censprob[i][l] = 0.0;
	}

	for (i = 0; i < m; i++)
		if ( GENOTYPE[M->J_P[i]][0]) // OBSERVED patterns
		{
			censprob[i][i] = Prob[i];
			for (l = 0; l < m; l++)
			{
				if (is_strict_after(M,l,i))
				{
					for (c=0; c < M->N_pa[l]; c++)
					{
						k = M->pa[l][c];  // index of c-th parent of i
						if(is_after(M,k,i))
						{
							j = M->pa_diff[l][c];  // index of differing event
							censprob[i][l] += (lambda[j] / lambda_exit[k]) * censprob[i][k];
						}
					}
				}
			}
		}
}


//Compute E[T|X] for all genotypes X
void compute_all_exp(double* lambda, model* M, int* g, double* lambda_exit, double* Prob, double** Exp, double** censprob, double*** censexp, int* lattice_index)
{
	/*
      Compute  E[X_i - max_{j \in \pa(i)} X_j | X compatible with Q]
      ( in code:  pos = i,  i = index over lattice elements (S) )
	 */

	int m = M->m;
	int n = M->n;

	int i, j, k, c, pred, pos, l;
	int all_pred_in_k;

	/*
  int* lattice_index = get_int_array(pow2(n+1));
  for (i=0; i < m; i++)
    lattice_index[M->J_P[i]]=i;
	 */

	for (pos=0;pos<=n;pos++)
		for(i=0;i<m;i++)
		{
			Exp[pos][i] = 0.0;
			for(l=0;l<m;l++)
				censexp[pos][i][l] = 0.0; // Clear everything
		}

#pragma omp parallel for private(i,c,k,j,g, all_pred_in_k, pred,l)
	for (pos=0; pos<=M->n; pos++)
	{
		for (i=0; i<M->m; i++) // First calculate events BEFORE censoring
		{
			Exp[pos][i] = 0.0;
			for (c=0; c<M->N_pa[i]; c++)
			{
				k = M->pa[i][c];  // index of c-th parent (in J_P) of i
				if (!GENOTYPE[M->J_P[k]][0] ) // BEFORE
				{
					j = M->pa_diff[i][c];  // index of differing event
					g = GENOTYPE[M->J_P[k]];

					all_pred_in_k = 1; // all predecessor events in J_P[k]?
					pred = 1;
					while (all_pred_in_k && (pred <= M->n))
					{
						all_pred_in_k = ((! M->P[pred][pos]) || (g[pred]));
						pred++;
					}

					Exp[pos][i] += (lambda[j] / lambda_exit[k]) * Exp[pos][k];
					if (all_pred_in_k)
					{
						Exp[pos][i] += ((1 - g[pos]) * lambda[j] * Prob[k]) / (lambda_exit[k] * lambda_exit[k]);
					}
				}
			}
		}

		for (i = 0; i < m; i++) // Now calculate the events AFTER censoring
			if ( GENOTYPE[M->J_P[i]][0] )
			{
				censexp[pos][i][i] = Exp[pos][i];
				for (l = 0; l < m; l++)
					if (is_strict_after(M,l,i))
					{
						for (c=0; c<M->N_pa[l]; c++)
						{
							k = M->pa[l][c];  // index of c-th parent (in J_P) of i
							if (is_after(M,k,i))
							{
								j = M->pa_diff[l][c];  // index of differing event
								g = GENOTYPE[M->J_P[k]];

								all_pred_in_k = 1; // all predecessor events in J_P[k]?
								pred = 1;
								while (all_pred_in_k && (pred <= M->n))
								{
									all_pred_in_k = ((! M->P[pred][pos]) || (g[pred]));
									pred++;
								}

								censexp[pos][i][l] += (lambda[j] / lambda_exit[k]) * censexp[pos][i][k];
								if (all_pred_in_k)
								{
									censexp[pos][i][l] += ((1 - g[pos]) * lambda[j] * censprob[i][k]) / (lambda_exit[k] * lambda_exit[k]);
								}
							}
						}
					}
				// Update Exp[pos][i] for s=1 (observation)
				Exp[pos][i] = censexp[pos][i][m-1];
			}
	}
}

// Compute the conditional probability Prob[X|Y]
void compute_condprob(model* M, data* D, int N_u, double* Prob, int only_falsepos, double epsilon, double** condprob)
{
	int n = M->n;
	int m = M->m;

	int i,k,d,index;

	double prob_tmp;

	for (k=0; k<N_u; k++)   // All patients
	{
		prob_tmp = 0;
		index = index_of(D[k].g, n+1);

		for (i=1; i<m; i++) // All patterns
		{
			condprob[i][k] = 0.0;
			// All observable genotypes or only false positives?
			//	if( (!only_falsepos && GENOTYPE[M->J_P[i]][0]) || (only_falsepos && M->J_P[i] >= index) )
			if( GENOTYPE[M->J_P[i]][0]) // Only OBSERVED patterns!
			{
				//d = ham_distance(M->J_P[i], index, &j, n+1);
				d = hamdist(M->J_P[i], index);
				condprob[i][k] = Prob[i] * power(epsilon, d) * power(1 - epsilon, n - d);
				prob_tmp += condprob[i][k];
			}
		}
		// Now divide by Prob[Y_k]
		for (i=1; i<m; i++) // All patterns
			if( GENOTYPE[M->J_P[i]][0])
				condprob[i][k] /= prob_tmp;
	}
}

// Compute E[T|Y], ie expected times T given observations Y
void compute_condexp(model* M, data* D, int N_u, double* Prob, double** condprob, double** Exp)
{
	int i,k,pos;
	int m = M->m;
	int n = M->n;

	for  (pos = 0; pos <= n; pos++) // Waiting times
		for (k = 0; k < N_u; k++) // Unique observations
		{
			D[k].t[pos] = 0.0;
			for (i = 0; i < m; i++) // Genotypes
				if( GENOTYPE[M->J_P[i]][0]) // Only OBSERVED patterns!
					D[k].t[pos] += Exp[pos][i] / Prob[i] * condprob[i][k];
		}
}

double EM_epsilon(model* M, data* D, int N_u, double* Prob, double** condprob, double* epsilon)
{
	int i, k, d, index;
	int iter = 0;
	int N = 0;
	int m = M->m;
	int n = M->n;

	double loglik = 0;
	double loglik_new = 0;
	double prob_tmp = 0;
	double H = 0;
	double delta_epsilon = *epsilon;

	// Calculate E[loglik(X,Y)|Y]
	while(iter < 10 && (iter < 2 || fabs(delta_epsilon)/ *epsilon > 1e-4 ) )
	{

		if ((iter > 5) && (loglik_new < loglik))  // should never happen!
		{
			fprintf(stderr, "Error in EM_epsilon: Likelihood is decreasing!\n");
			exit(1);
		}

		loglik = loglik_new;
		loglik_new = 0;
		prob_tmp =0;

		/* E-step, compute expected Hamming distance E[H(X,Y)|Y] */
		compute_condprob(M, D, N_u, Prob, 0, *epsilon, condprob);

		H = 0;
		N = 0;
		for (k = 0; k < N_u; k++)
		{
			index = index_of(D[k].g, n+1);
			N += D[k].count;
			for(i = 1; i < m; i++)
				if( GENOTYPE[M->J_P[i]][0])
				{
					d = hamdist(M->J_P[i], index);
					H += d * condprob[i][k] * D[k].count;
					prob_tmp += Prob[i] * power(*epsilon, d) * power(1 - *epsilon, n - d);
				}
			loglik_new += log (prob_tmp) * D[k].count;
			prob_tmp = 0.0;
		}
		H /= N;

		/* M-step */
		delta_epsilon = H/n - *epsilon;
		*epsilon = H / n;
		iter++;
	}
	return loglik_new;
}

int compute_loglik(model* M, data* D, int N_u, double* Prob, double** condprob, double* epsilon, double* loglik){
	int i, k, d, index;
	int m = M->m;
	int n = M->n;

	double prob_tmp = 0;
	compute_condprob(M, D, N_u, Prob, 0, *epsilon, condprob);

	for (k = 0; k < N_u; k++)
	{
		index = index_of(D[k].g, n+1);
		for(i = 1; i < m; i++)
			if( GENOTYPE[M->J_P[i]][0])
			{
				d = hamdist(M->J_P[i], index);
				prob_tmp += Prob[i] * power(*epsilon, d) * power(1 - *epsilon, n - d);
			}
		loglik[k] = log (prob_tmp) ;
		prob_tmp = 0.0;
	}
	return 0;
}

double compute_total_loglik(model* M, data* D, int N_u, double* Prob, double** condprob, double* epsilon){
	int k;
	double total_loglik = 0.0;
	double* loglik = get_double_array(N_u);
	compute_loglik(M, D, N_u, Prob, condprob, epsilon, loglik);
	for (k = 0; k < N_u; k++)
		total_loglik += loglik[k] * D[k].count;
	return total_loglik;
}

void compute_lambda_exit(double *lambda, model* M, double* lambda_exit)
{
	int i, j, c;

	for (i=0; i<M->m; i++)
	{
		lambda_exit[i] = 0.0;

		for (c=0; c<M->N_ch[i]; c++)
		{
			j = M->ch_diff[i][c];  // index of differing event
			lambda_exit[i] += lambda[j];

		}
	}

}

double compute_GPS ( double* lambda, model* M, int* sublattice, double* lambda_exit, double* Prob, double* Exp )
{
	/*
      Compute  E[max{X_i}_{i\in S}  | X_i < X_j \forall j \in P\S]
      ( i = index over sublattice elements (S), S = sublattice )
	 */
	double GPS = 0;
	int i, j, k, c;

	for (i=0; i<M->m; i++)
	{
		Exp[i] = 0.0;

		if (sublattice[i] && !GENOTYPE[M->J_P[i]][0])
		{
			for (c=0; c<M->N_pa[i]; c++)
			{
				k = M->pa[i][c];  // index of c-th parent (in J_P) of i
				if (sublattice[k] && !GENOTYPE[M->J_P[k]][0])
				{
					j = M->pa_diff[i][c];  // index of differing event

					Exp[i] += (lambda[j] / lambda_exit[k]) * Exp[k] + (lambda[j] / (lambda_exit[k] * lambda_exit[k])) * Prob[k];

					GPS = Exp[i] / Prob[i];
				}
			}
		}
	}
	return GPS;
}

void compute_all_GPS ( double* lambda, model* M, double* lambda_exit, double* Prob, int* lattice_index, double* all_GPS )
{
	/*Compute essentially the same as compute_GPS() yet for all possible observations*/

	int i, j, k, c;
	double* Exp = get_double_array(M->m);
	int power = pow2(M->n);

	for (i=0; i<M->m; i++)
	{
		Exp[i] = 0.0;

		if ( !GENOTYPE[M->J_P[i]][0])
		{
			for (c=0; c<M->N_pa[i]; c++)
			{
				k = M->pa[i][c];  // index of c-th parent (in J_P) of i
				if ( !GENOTYPE[M->J_P[k]][0])
				{
					j = M->pa_diff[i][c];  // index of differing event

					Exp[i] += (lambda[j] / lambda_exit[k]) * Exp[k] + (lambda[j] / (lambda_exit[k] * lambda_exit[k])) * Prob[k];

					all_GPS[i] = Exp[i] / Prob[i];
				}
			}
			Exp[lattice_index[M->J_P[i]+power]] = (lambda[0] / lambda_exit[i]) * Exp[i] + (lambda[0] / (lambda_exit[i] * lambda_exit[i])) * Prob[i];
			all_GPS[lattice_index[M->J_P[i]+power]] = Exp[lattice_index[M->J_P[i]+power]] / Prob[lattice_index[M->J_P[i]+power]];
		}
	}
	free(Exp);
}

void compute_cond_GPS (model* M, data* D, int N_u, double* Prob, double** condprob, double* all_GPS, double* cond_GPS)
{
	int i,k;

	for (k = 0; k < N_u; k++) // Unique observations
	{
		cond_GPS[k] = 0.0;
		for (i = 0; i < M->m; i++) // Genotypes
			if( GENOTYPE[M->J_P[i]][0]) // Only OBSERVED patterns!
				cond_GPS[k] += all_GPS[i] * condprob[i][k];
	}
}

void GPS (model* M, data*D, int N_u, double* lambda, double epsilon, double* all_GPS, double* cond_GPS)
{
	int i;
	int m = M->m;
	int n = M->n;

	double* lambda_exit = get_double_array(m);
	double* Prob = get_double_array(m);
	double** condprob = get_double_matrix(m, N_u);

	int* lattice_index = get_int_array(pow2(n+1));
	for (i=0; i < m; i++)
		lattice_index[M->J_P[i]]=i;

	compute_lambda_exit(lambda, M, lambda_exit);
	compute_all_prob(lambda, M, lambda_exit, Prob, lattice_index);
	compute_condprob(M, D, N_u, Prob, 0, epsilon, condprob);

	// Compute GPS for all observable genotypes
	compute_all_GPS (lambda, M, lambda_exit, Prob, lattice_index, all_GPS);

	// Compute GPS given the observations
	compute_cond_GPS(M, D, N_u, Prob, condprob, all_GPS, cond_GPS);

	free(Prob);
	free(lambda_exit);

	for(i=0;i<m;i++)
	{
		free(condprob[i]);
	}
	free(condprob);
	free(lattice_index);
	free(all_GPS);
}

/* Compute the most likely hidden state X*/
void compute_hidden_patterns(model* M, data*D, int N_u, double* lambda, double epsilon, double** condprob, int* ml_pat, double** exp_pat)
{
	int i,j,k;
	for (k=0; k<N_u; k++)
	{
		ml_pat[k] = 0;
		for (i=0; i<M->m; i++)
		{
			// Find MAP genotype \tilde{X}
			if (condprob[i][k] > condprob[ml_pat[k]][k])
				ml_pat[k] = i;
			// Compute E[X|Y_k]
			for (j = 0; j < M->n; j++)
				exp_pat[j][k] += GENOTYPE[M->J_P[i]][j + 1] * condprob[i][k];
		}
	}
}

/* 14-Jul-2009 Compute Prob[Y] for all Y in {0,1}^n */
double* compute_ProbY(model* M, double* Prob, int only_falsepos, double epsilon)
{
	int n = M->n;
	int m = M->m;

	int i,k,d;

	int size = pow2(n);
	double* ProbY = get_double_array(size);

	for (k=0; k < size; k++)   // All possible observations Y
	{
		ProbY[k] = 0.0;
		for (i=0; i<m; i++) // All patterns
		{
			// All observable genotypes or only false positives?
			//	if( (!only_falsepos && GENOTYPE[M->J_P[i]][0]) || (only_falsepos && M->J_P[i] >= index) )
			if( GENOTYPE[M->J_P[i]][0]) // Only OBSERVED patterns!
			{
				d = hamdist(M->J_P[i], k+pow2(n));
				ProbY[k] += Prob[i] * power(epsilon, d) * power(1 - epsilon, n - d);
			}
		}
	}
	return ProbY;
}

/* 14-Jul-2009 Compute Prob[Y] in ct-cbn model */
double* compute_ProbY_ctcbn(model* M, double* Prob)
{
	int n = M->n;
	int m = M->m;

	int i,k;

	int size = pow2(n);
	double* ProbY_ctcbn = get_double_array(size);
	int* is_compatible = get_int_array(size);
	int* which = get_int_array(size);
	int N_compatible = 0;

	for (k=0; k < size; k++){ // All possible observations Y
		i = 0;
		while (!is_compatible[k] && i < m){ // All patterns
			if ( M->J_P[i] == k + pow2(n)){
				is_compatible[k] = 1;
				which[k] = i;
				N_compatible++;
			}
			i++;
		}
	}

	double alpha = (double) N_compatible / (double) size;
	double q_eps = 1.0 / (size - N_compatible);  // "+1" because n is the number of events,


	for (k=0; k < size; k++){   // All possible observations Y
		if ( is_compatible[k] )
			ProbY_ctcbn[k] = alpha * Prob[which[k]];
		else
			ProbY_ctcbn[k] = (1 - alpha) * q_eps;
	}

	return ProbY_ctcbn;
}

void poset_from_data(int** P, int n, int* pat, int** Q)
{
	int i, j;

	for (i=0; i<=n; i++)
		for (j=0; j<=n; j++)
			Q[i][j] = P[i][j];

	for (i=1; i<=n; i++)
	{
		Q[i][0] = pat[i];
		Q[0][i] = (1 + pat[i]) % 2;
	}

}



void refine_poset(int** P, const int n, double* t, int** Q)
{
	int i, j;

	for (i=0; i<=n; i++)
		for (j=0; j<=n; j++)
			Q[i][j] = P[i][j];

	for (i=1; i<=n; i++)
	{
		if (t[i] < t[0])
			Q[i][0] = 1;
		else
			Q[0][i] = 1;
	}

}



void sublattice(int* J_P, const int m, int** Q, const int n, int* J_Q)
{
	int i, j, l;
	int *g;
	int brk;

	for (l=0; l<m; l++)
	{
		brk = 0;
		g = GENOTYPE[J_P[l]];

		// check compatibility:
		J_Q[l] = 1;
		for (i=0; i<n && !brk; i++)
			for (j=0; j<n; j++)
			{
				if (Q[i][j] && (! g[i]) && g[j])
					/*  i.e.,  i < j  in Q,  but  g_i==0  and  g_j==1  */
				{
					J_Q[l] = 0;  // incompatible!
					brk = 1;
					break;
				}
			}
	}

}



double observed_log_likelihood(int** P, int n, double** t, int N, double* lambda)
{
	// assumes cover relations in P

	int i, j, k;
	double max, loglik = 0.0;

	double* log_lambda = get_double_array(n+1);
	for (j=0; j<=n; j++)
		log_lambda[j] = log(lambda[j]);

	for (i=1; i<=n; i++)
	{
		for (k=0; k<N; k++)
		{
			max = 0.0;
			for (j=1; j<=n; j++)
			{
				if (P[j][i])  // i.e.  j covers i
					max = MAX(max, t[k][j]);
			}
			loglik += log_lambda[i] - lambda[i] * (t[k][i] - max);
		}
	}
	free(log_lambda);

	return loglik;
}



double hidden_log_likelihood(int** P, int n, int* J_P, int m, int** pat, int N, double* lambda)
{
	double loglik = 0.0;

	// not implemented, see EM!

	return loglik;
}



void MLE(int n, data* D, int N_u, double* lambda)
{
	/* Maximum likelihood estimation of lambda */

	int i, k, N;
	double sum;

	for (i=0; i<=n; i++)
	{
		sum = 0.0;
		N = 0;
		for (k=0; k<N_u; k++)
			if (D[k].is_compatible)
			{
				N += D[k].count;
				sum += D[k].count * D[k].t[i];
			}
		lambda[i] = (double) N / sum;
	}

}

void MLE2(int n, data* D, int N_u, double* lambda)
{
	/* Maximum likelihood estimation of lambda */

	int i, k, N;
	double sum;

	for (i=0; i<=n; i++)
	{
		sum = 0.0;
		N = 0;
		for (k=0; k<N_u; k++)
			//	if (D[k].is_compatible)
		{
			N += D[k].count;
			sum += D[k].count * D[k].t[i];
		}
		lambda[i] = (double) N / sum;
	}

}


double norm_diff(double* x, double* y, int n)
{
	int j;
	double summand, sum = 0.0;

	for (j=0; j<n; j++)
	{
		summand = (x[j] - y[j]) / x[j];
		sum += summand * summand;
	}

	return sum;
}



double EM(model* M, data* D, int N_u, double* lambda, double* alpha, double q_eps)
{
	int j, k;
	int iter = 0;
	double S = lambda[0];
	double loglik = 0.0;
	double loglik_new = 0.0;
	double log_alpha = log(*alpha);

	int m = M->m;
	int n = M->n;

	int* g = get_int_array(n+1);
	double* lambda_new = get_double_array(n+1);
	double* lambda_exit = get_double_array(m);
	double* Prob = get_double_array(m);
	double** Exp = get_double_matrix(n+1, m);

	while ((norm_diff(lambda, lambda_new, n+1) > EM_ACCURACY || (iter < 2)) && (iter < EM_MAX_ITER))
	{
		if ((iter > 2) && (loglik_new < loglik))  // should never happen!
		{
			fprintf(stderr, "Error: Likelihood is decreasing!\n");
			exit(1);
		}

		loglik = loglik_new;
		loglik_new = 0.0;
		for (j=0; j<=n; j++)
			lambda_new[j] = lambda[j];

		/*****************************************************************
	    E step:  compute  E[X_i - max_{j \in pa(i)} X_j | lambda, Q]
		 ******************************************************************/

		compute_lambda_exit(lambda, M, lambda_exit);

		for (k=0; k<N_u; k++)
			if (D[k].is_compatible)
			{
				compute_prob(lambda, M, D[k].J_Q, lambda_exit, Prob);

				if (Prob[m-1] <= 0.0)  // should never happen!
				{
					fprintf(stderr, "Error: Genotype %d incompatible with poset!\n", k);
					exit(1);
				}

				loglik_new += (D[k].count * (log_alpha + log(Prob[m-1])));

				compute_exp(lambda, M, g, D[k].J_Q, lambda_exit, Prob, Exp);

				for (j=0; j<=n; j++)
					D[k].t[j] = Exp[j][m-1] / Prob[m-1];
			}

		if (verbose)
		{
			printf("%d\t", iter);
			printf(DOUBLE_FORMAT, loglik_new);  printf("\t");
			print_double_array(lambda, n+1);
		}

		/************************************
	    M step:  compute  lambda_ML(t)
		 *************************************/

		MLE(n, D, N_u, lambda);

		double fac = S / lambda[0];
		for (j=0; j<=n; j++)
			lambda[j] *= fac;

		iter++;
	}

	if (iter >= EM_MAX_ITER)
		fprintf(stderr, "Warning: EM aborted without convergence after %d iterations!\n", iter);


	if (write_expected_times) {  // for diagnostic purposes only

		char filestem[15] = "X";
		char suffix[15] = ".tab";
		char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
		strcat(filename, filestem);
		strcat(filename, suffix);
		FILE *output;
		if ( (output = fopen(filename, "w")) == NULL )
		{
			fprintf(stderr, "Error:  Could not write to file %s\n", filename);
			exit(1);
		}

		fprintf(output, "count\tis_comp\tresp");
		for (j=0; j<=n; j++)
			fprintf(output, "\tX_%d", j);
		for (j=0; j<=n; j++)
			fprintf(output, "\tt_%d", j);
		fprintf(output, "\n");


		for (k=0; k<N_u; k++)
		{
			fprintf(output, "%d", D[k].count);
			fprintf(output, "\t%d", D[k].is_compatible);

			double resp = 0.0; // responsibility of poset for observation
			if (*alpha == 1.0)
				resp = 1.0;
			else
			{
				compute_prob(lambda, M, D[k].J_Q, lambda_exit, Prob);
				if (Prob[m-1] > 0.0)
					resp = (*alpha * Prob[m-1]) / ((1.0 - *alpha)*q_eps + *alpha*Prob[m-1]);
			}

			fprintf(output, "\t%g", resp);

			for (j=0; j<=n; j++)
				fprintf(output, "\t%d", D[k].g[j]);
			for (j=0; j<=n; j++)
				fprintf(output, "\t%g", D[k].t[j]);
			fprintf(output, "\n");

		}

		fclose(output);

	}


	free(lambda_new);
	for (j=0; j<=n; j++)
		free(Exp[j]);
	free(Exp);
	free(Prob);
	free(lambda_exit);
	free(g);

	return loglik;
}


// CALCULATE NESTED EM
double EM_EM(model* M, data* D, int N_u, double* lambda, double* epsilon)
{
	//int i, k,l;
	int i,j;
	int iter = 0;
	double S = lambda[0];
	double loglik = 0.0;
	double loglik_new = 0.0;
	//double log_alpha = log(*alpha);

	int m = M->m;
	int n = M->n;

	int* g = get_int_array(n+1);
	double* lambda_new = get_double_array(n+1);
	double* lambda_exit = get_double_array(m);
	double* Prob = get_double_array(m);
	double** Exp = get_double_matrix(n+1, m);
	double** condprob = get_double_matrix(m, N_u);
	double** censprob = get_double_matrix(m,m);
	double*** censexp = get_double_cube(n+1,m,m);

	int* lattice_index = get_int_array(pow2(n+1));
	for (i=0; i < m; i++)
		lattice_index[M->J_P[i]]=i;

	// Estimate LAMBDA
	while ((norm_diff(lambda, lambda_new, n+1) > EM_ACCURACY || (iter < 2)) && (iter < EM_MAX_ITER) && (loglik_new - loglik > 1e-4 || (iter < 2)))
	{

		if ((iter > 5) && (loglik_new < loglik))  // should never happen!
		{
			fprintf(stderr, "Error in EM_EM: Likelihood is decreasing!\n");
			exit(1);
		}


		loglik = loglik_new;
		loglik_new = 0.0;
		for (j=0; j<=n; j++)
			lambda_new[j] = lambda[j];

		/*****************************************************************
	    E step:  compute  E[X_i - max_{j \in pa(i)} X_j | lambda, Q]
		 ******************************************************************/
		compute_lambda_exit(lambda, M, lambda_exit);

		// Calculate all probabilities given lambda_i
		compute_all_prob(lambda, M, lambda_exit, Prob, lattice_index);

		/*
	if (Prob[m-1] <= 0.0)  // should never happen!
	{
	fprintf(stderr, "Error: Genotype %d incompatible with poset!\n", k);
	exit(1);
	}
		 */

		// Estimate epsilon, return loglik
		loglik_new = EM_epsilon(M, D, N_u, Prob, condprob, epsilon);

		/* CONDPROB OK.
      for (i=0;i<m;i++)
	if(GENOTYPE[M->J_P[i]][0])
	  {
	    for(k=0;k<N_u;k++)
	      printf("%f\t",condprob[i][k]);
	    printf("\n");
	  }
      printf("\n");
		 */

		compute_censored_prob(lambda, M, lambda_exit, Prob, censprob, lattice_index);

		/* SEEMS OK...
      for(i=0; i<m;i++)
	if(GENOTYPE[M->J_P[i]][0])
	  {
	    print_int_array(GENOTYPE[M->J_P[i]],n+1);
	    printf("%f\n", Prob[i]);
	    printf("---\n");
	    for(l = 0; l < m; l++)
	      if( M->J_P[l] >= M->J_P[i] &&  GENOTYPE[M->J_P[l]][0])
		{
		  print_int_array(GENOTYPE[M->J_P[l]],n+1);
		  printf("%f\n",censprob[i][l]);
		}
	    printf("\n");
	  }
      printf("\n");
		 */

		// Compute all exp waiting times E[T|S]
		compute_all_exp(lambda, M, g, lambda_exit, Prob, Exp, censprob, censexp, lattice_index);

		// Compute expected waiting times E[T|Y]
		compute_condexp( M, D, N_u, Prob, condprob, Exp);

		if (verbose)
		{
			printf("%d\t%f\t", iter,*epsilon);
			printf(DOUBLE_FORMAT, loglik_new);  printf("\t");
			print_double_array(lambda, n+1);
		}
		fflush(stdout);

		/************************************
	    M step:  compute  lambda_ML(t)
		 *************************************/

		MLE2(n, D, N_u, lambda);

		double fac = S / lambda[0];
		for (j=0; j<=n; j++)
			lambda[j] *= fac;

		iter++;
	}

	if (iter >= EM_MAX_ITER)
		fprintf(stderr, "Warning: EM aborted without convergence after %d iterations!\n", iter);

	/*
  if (write_expected_times) {  // for diagnostic purposes only

    char filestem[15] = "X";
    char suffix[15] = ".tab";
    char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
    strcat(filename, filestem);
    strcat(filename, suffix);
    FILE *output;
    if ( (output = fopen(filename, "w")) == NULL )
      {
	fprintf(stderr, "Error:  Could not write to file %s\n", filename);
	exit(1);
      }

    fprintf(output, "count\tis_comp\tresp");
    for (j=0; j<=n; j++)
      fprintf(output, "\tX_%d", j);
    for (j=0; j<=n; j++)
      fprintf(output, "\tt_%d", j);
    fprintf(output, "\n");


    for (k=0; k<N_u; k++)
      {
	fprintf(output, "%d", D[k].count);
	fprintf(output, "\t%d", D[k].is_compatible);

	double resp = 0.0; // responsibility of poset for observation
	if (*alpha == 1.0)
	  resp = 1.0;
	else
	  {
	    compute_prob(lambda, M, D[k].J_Q, lambda_exit, Prob);
	    if (Prob[m-1] > 0.0)
	      resp = (*alpha * Prob[m-1]) / ((1.0 - *alpha)*q_eps + *alpha*Prob[m-1]);
	  }

	fprintf(output, "\t%g", resp);

	for (j=0; j<=n; j++)
	  fprintf(output, "\t%d", D[k].g[j]);
	for (j=0; j<=n; j++)
	  fprintf(output, "\t%g", D[k].t[j]);
	fprintf(output, "\n");

      }

    fclose(output);


  }
	 */



	free(lambda_new);
	for (j=0; j<=n; j++)
	{
		free(Exp[j]);
		for(i=0;i<m;i++)
			free(censexp[j][i]);
		free(censexp[j]);
	}
	free(censexp);
	free(Exp);

	free(Prob);
	free(lambda_exit);
	free(g);

	for(i=0;i<m;i++)
	{
		free(censprob[i]);
		free(condprob[i]);
	}
	free(censprob);
	free(condprob);
	free(lattice_index);
	return loglik;
}

void guess_lambda(int** P, int n, data* D, int N_u, double* lambda)
{
	int i, j, k;
	double S = lambda[0];

	for (j=1; j<=n; j++)
	{
		int total_count = 2;  // pseudo count
		int mut_count = 1;  // pseudo count
		/* Pseudo counts should ensure that total_counts > 0 and p < 1 (see below).
	 Note that p==1 means that the event is redundant. */
		for (k=0; k<N_u; k++)
		{
			int parents_ok = 1;
			for (i=1; i<n; i++)
			{
				if (P[i][j] && (D[k].g[i] == 0))
				{  // ie. i covers j, but sample k lacks mutation i
					parents_ok = 0;
					break;
				}
			}
			if (parents_ok)
			{
				total_count += D[k].count;
				mut_count += D[k].count * D[k].g[j];
			}
		}
		double p = (double) mut_count / (double) total_count;
		lambda[j] = (p / (1.0 - p)) * S;
	}

}



double estimate_model_parameters(model* M, data* D, int N_u, int R, double* lambda_opt, double* alpha)
{
	int j, r, k;
	double loglik = 0.0;
	double loglik_opt = 0.0;
	double* lambda = get_double_array(M->n+1);

	// compute alpha_ML and log-likelihood contribution from incompatible samples
	int N = 0;
	int N_compatible = 0;
	for (k=0; k<N_u; k++)
	{
		N += D[k].count;
		N_compatible += D[k].is_compatible * D[k].count;
	}
	*alpha = (double) N_compatible / (double) N;
	double q_eps = 1.0 / (pow2(M->n + 1) - M->m);  // "+1" because n is the number of events,
	// but P_eps has n+1 elements, namely [n] U {s}.
	double loglike_incompatible = 0.0;
	if (N > N_compatible)
		loglike_incompatible = (N - N_compatible) * (log(1.0 - *alpha) + log(q_eps));

	// compute lambda_ML by running the EM
	lambda[0] = lambda_opt[0];
	for (r=0; r<R; r++)  // EM repetitions
	{
		guess_lambda(M->P, M->n, D, N_u, lambda);
		loglik = EM(M, D, N_u, lambda, alpha, q_eps) + loglike_incompatible;
		if (verbose)  printf("\n");
		if ((r == 0) || (loglik > loglik_opt))
		{
			loglik_opt = loglik;
			for (j=1; j<=M->n; j++)
				lambda_opt[j] = lambda[j];
		}
	}

	free(lambda);

	return loglik_opt;
}



void draw_samples(int** P, double* lambda, int* lin_ext, int n, int** pat, double** t, int N)
{
	int i, j, k, idx;
	double max;
	double* Z = get_double_array(n+1);

	for (k=0; k<N; k++)
	{
		for (j=0; j<=n; j++)
			Z[j] = gsl_ran_exponential(RNG, 1/lambda[j]);

		t[k][0] = Z[0];
		for (i=0; i<n; i++)
		{
			idx = lin_ext[i];  // visit events in order of a linear extension of the poset
			max = 0.0;
			for (j=1; j<=n; j++)
			{
				if (P[j][idx])  // j covers idx
					// Note: j has been visited before because of linear extension order
					max = MAX(max, t[k][j]);
			}
			t[k][idx] = max + Z[idx];
		}

		pat[k][0] = 1;
		for (j=1; j<=n; j++)
			if (t[k][j] > t[k][0])
				pat[k][j] = 0;  // not observed
			else
				pat[k][j] = 1;  // observed
	}

	free(Z);
}



data* make_data_set(int** pat, int N, int n, int* N_u, int* pat_idx)
{
	int j, k, l, skip;

	int* idx = get_int_array(N);
	int* count = get_int_array(N);

	for (k=0; k<N; k++)
		idx[k] = index_of(pat[k], n+1);

	// count unique patterns:
	int c = 0;
	for (k=0; k<N; k++)
	{
		skip = 0;
		for (l=0; l<k; l++)
		{
			if (idx[k] == idx[l]) // observed before?
			{
				count[l]++;
				pat_idx[k] = pat_idx[l];
				skip = 1;
				break;
			}
		}
		if (! skip) // Unique!
		{
			count[k]++;
			pat_idx[k] = c;
			c++;
		}
	}

	*N_u = 0;
	for (k=0; k<N; k++)
		*N_u += (count[k] > 0);

	if (verbose)  printf("N_u = %d unique patterns\n", *N_u);

	data* D = calloc(*N_u, sizeof(data));

	c = 0;
	for (k=0; k<N; k++)
	{
		//    pat_idx[k] = c;
		if (count[k] > 0)
		{
			//	  pat_idx[k]=c;
			D[c].count = count[k];
			D[c].g = get_int_array(n+1);
			for (j=0; j<=n; j++)
				D[c].g[j] = pat[k][j];
			c++;
		}
	}

	return D;
}



void construct_sublattices(data* D, int N_u, model* M)
{
	int k;
	int m = M->m;
	int n = M->n;

	for (k=0; k<N_u; k++)
	{
		D[k].Q = get_int_matrix(n+1, n+1);
		poset_from_data(M->P, n, D[k].g, D[k].Q);

		D[k].J_Q = get_int_array(m);
		sublattice(M->J_P, m, D[k].Q, n+1, D[k].J_Q);

		D[k].t = get_double_array(n+1);  // times: missing data
	}

}



int lt_poset(int a, int b, int**P, int n)
{
	// check if a < b in P

	int i, j;
	int* added = get_int_array(n+1);

	queue q;
	init_queue(&q);
	enqueue(&q, b);
	added[b] = 1;

	// visit below(b) and check for a:
	while (empty(&q) == FALSE) {
		j = dequeue(&q);
		for (i=0; i<=n; i++)
			if (P[i][j] && (! added[i]))  // i.e. i < j
			{
				if (i == a)  // have it!
				{
					free(added);
					return 1;
				}
				enqueue(&q, i);
				added[i] = 1;
			}
	}

	free(added);
	return 0;
}



int compare_violation_pairs (const void *A, const void *B)
{
	int* a = *(int**) A;
	int* b = *(int**) B;

	// random tie breaking:
	double da = (double) a[2] + (double) rand() / (double) RAND_MAX;
	double db = (double) b[2] + (double) rand() / (double) RAND_MAX;

	return (da - db);
}



int** violation_map(data* D, int N_u, int n)
{
	int i, j, k, idx;

	int** V = get_int_matrix((n+1) * (n+1), 3);

	// count violations
	idx = 0;
	for (i=1; i<=n; i++)
		for (j=1; j<=n; j++)
		{
			V[idx][0] = i;
			V[idx][1] = j;
			V[idx][2] = 0;  // counts how often i<j is violated by the data
			for (k=0; k<N_u; k++)
			{
				if ((! D[k].g[i]) && D[k].g[j])  // i.e., i < j but g_i==0 and g_j==1
					V[idx][2] += D[k].count;
			}
			idx++;
		}

	qsort(V, idx, sizeof(int *), compare_violation_pairs);  // small violators first

	/*  This is how the comparison function is called:
      int c05 = compare_violation_pairs(&V[0], &V[5]);  */

	return V;
}

// Changed by Moritz 09/04/2008
int reduce_to_cover_relations ( int** P, int n )
{
	int i,j,k;
	queue q;
	int stat = 0;

	// Sample all nodes
	for (i=1; i<=n; i++)
	{
		int* visit = get_int_array(n+1);
		init_queue(&q);

		// Fill queue with children of i
		for (j=1; j<=n; j++)
			if (P[i][j])
			{
				enqueue(&q,j);
			}

		// Walk through grandchildren
		while (empty(&q) == FALSE)
		{
			j = dequeue(&q);
			for (k=1; k<=n; k++)
			{
				if (P[j][k] && !(visit[k]))
				{
					visit[k] = 1;
					enqueue(&q,k);

					// Remove non-cover relations
					if (P[i][k])
					{
						P[i][k] = 0;
						stat = 1; // Report changes
					}

					// Check if cyclic
					if (P[k][i])
						return 2; // Fatal
				}
			}
		}
		free(visit);
	}
	return stat;
}

void maximal_poset(data* D, int N_u, int n, double eps, int** P)
{
	/*
      Construct maximal poset by including all relations
      that violate eps*100% or less of the data
	 */

	int i,j , k, idx;
	int N_pairs = (n+1) * (n+1);

	int N = 0;
	for (k=0; k<N_u; k++)
		N += D[k].count;

	int** V = violation_map(D, N_u, n);

	// build poset:
	for (i=0; i<=n; i++)
		for (j=0; j<=n; j++)
			P[i][j] = 0;

	for(idx=0; idx<N_pairs; idx++)
	{
		i = V[idx][0];
		j = V[idx][1];
		if (i && j)
			if ((double) V[idx][2] / (double) N <= eps)
				if (! lt_poset(j, i, P, n))
					P[i][j] = 1;
	}

	for (i=0; i<N_pairs; i++)
		free(V[i]);
	free(V);

	// reduce to cover relations:
	for (i=1; i<=n; i++)
		P[i][i] = 0;

	reduce_to_cover_relations (P, n);

	/* Changed by Moritz 09/04/2008
  int changed = 1;
  while (changed)
    {
      changed = 0;
      for (i=1; i<=n; i++)
	for (j=1; j<=n; j++)
	  if (P[i][j])
	    for (k=1; k<=n; k++)
	      if (P[i][k] && P[k][j])
		{
		  P[i][j] = 0;
		  changed = 1;
		}
    }
	 */
}

void compatibility(data* D, int N_u, model* M)
{
	int i, j, k;

	for (k=0; k<N_u; k++)
	{
		D[k].is_compatible = 1;
		i = 1;
		while ((i <= M->n) && D[k].is_compatible)
		{
			j = 1;
			while ((j <= M->n) && D[k].is_compatible)
			{
				if (M->P[i][j] && (! D[k].g[i]) && D[k].g[j])
					/* It is sufficient to check the cover relations! */
					//if (lt_poset(i, j, M->P, M->n) && (! D[k].g[i]) && D[k].g[j])
					D[k].is_compatible = 0;
				j++;
			}
			i++;
		}
	}

}



void select_poset(int k, double eps, model* M, double* lambda, data* D, int N_u, int R, int mode, int print)
{

	double loglik = 0.0;
	double alpha;  // mixing parameter, fraction of compatible genotypes

	if (mode != LEARN_PARAM)
		maximal_poset(D, N_u, M->n, eps, M->P);

	compatibility(D, N_u, M);

	M->J_P = bfs_order_ideals(M->P, M->n+1, &(M->m), M->lin_ext);
	parents(M);
	children(M);
	if (verbose)  print_model(M);

	construct_sublattices(D, N_u, M);
	if (verbose)  print_data(D, N_u, M->n, M->m);

	if (mode != LEARN_POSET)
	{
		loglik = estimate_model_parameters(M, D, N_u, R, lambda, &alpha);
		if(print){
			printf("%d\t%g\t%g\t%g\t", k, eps, alpha, loglik);
			print_double_array(lambda, M->n+1);
		}
	}

	//  free_lattice(M);

}


/* LOCAL SEARCH ALGORITHM SECTION*/
void make_model(model* M, int n, int** P)
{
	int i,j;

	M->n = n;
	M->P = get_int_matrix(n+1, n+1);
	precompute_binary(M->n+1);
	M->lin_ext = get_int_array(M->n);  // a linear extension of the poset

	for (i=0;i<n+1;i++)
		for (j=0;j<n+1;j++)
			M->P[i][j] = P[i][j];

	M->J_P = bfs_order_ideals(M->P, M->n+1, &(M->m), M->lin_ext);
	parents(M);
	children(M);
}

double try_edge(model* M, model* M2, data* D, int N_u, double* lambda, double* epsilon, double loglik, double** loglik_next, double T)
{
	int i,j,k,c;
	int reject = 0;

	int n = M->n;

	int R1,R2;
	int* R4 = get_int_array(n*n);
	double* R5 = get_double_array(n*n);
	double R3;

	double loglik_new = loglik;
	double* lambda_new = get_double_array(n+1);
	lambda_new[0] = lambda[0];
	double boltz;

	int N_compatible = 0;
	int N = 0;
	double alpha, alpha_new;
	double** alpha_all = get_double_matrix(n,n);

	compatibility(D, N_u, M);
	for (k=0; k<N_u; k++)
	{
		N += D[k].count;
		N_compatible += D[k].is_compatible * D[k].count;
	}
	alpha = (double) N_compatible / (double) N;
	alpha_new = alpha;

	/* Draw n^2 random numbers */
	for (i=0;i<n*n;i++)
		R5[i] = (double) rand() / (double) RAND_MAX;

	/* Sort to generate integers */
	for (i=0;i<n*n;i++)
	{
		c =0;
		for(j=0;j<n*n;j++)
			if(R5[i]>R5[j])
				c++;
		R4[i] = c;
	}

	int iter = 0;
	int poset_stat;

	/* Compute no. of violated data */
	while (iter < n*n)
	{
		reject = 0;

		/* Copy poset*/
		for (i=1;i<n+1;i++)
			for(j=1;j<n+1;j++)
				M2->P[i][j] = M->P[i][j];

		R1 = (int) iter / n + 1;
		R2 = iter % n + 1;

		//      printf("%i, %i\n", R1, R2);

		if (R1 == R2)
		{
			reject = 1;
			alpha_all[R1-1][R2-1] = 0.0;
		}

		if (!reject)
		{
			M2->P[R1][R2] = 1 - M2->P[R1][R2]; // Change edge

			if (M2->P[R1][R2]) // added edge (needs testing)
			{
				//printf("Added.\n");
				poset_stat = reduce_to_cover_relations (M2->P, n);
				if (poset_stat != 2) // Fatal
				{
					if (!M2->P[R1][R2]) // ie, suggested relation has been removed
					{
						c = 0;
						for(i=1;i<n+1;i++)
							if(M2->P[R1][i] && M2->P[i][R2]) // Check if single intermediate state, ie R1 -> j -> R2
							{
								j = i;
								c++;
							}
						if(c == 1)
						{
							M2->P[j][R2] = 0; // Remove j -> R2
							M2->P[R1][R2] = 1; // Add R1 -> R2
							if (reduce_to_cover_relations (M2->P, n) > 0)
								reject = 1;
						}
						else reject = 1;
					}
					compatibility(D, N_u, M2);
					N_compatible = 0;
					for (k=0; k<N_u; k++)
					{
						N_compatible += D[k].is_compatible * D[k].count;
					}
					alpha_all[R1-1][R2-1] = (double) N_compatible / (double) N;
				}
				else
					reject = 1;
			}

			else // deleted edge
			{
				//printf("Deleted.\n");
				poset_stat = reduce_to_cover_relations (M2->P, n);
				compatibility(D, N_u, M2);
				N_compatible = 0;
				for (k=0; k<N_u; k++)
				{
					N_compatible += D[k].is_compatible * D[k].count;
				}
				alpha_all[R1-1][R2-1] = (double) N_compatible / (double) N;
			}
		}

		if (reject)
			alpha_all[R1-1][R2-1] = 0.0;

		iter++;
	}

	//print_double_matrix(alpha_all, n, n);
	//return;

	/*   int iter = 0; */
	/*   int poset_stat; */
	/*   while (iter < n*n) */
	/*     { */
	/*       reject = 0; */

	/*       /\* Copy poset*\/ */
	/*       for (i=1;i<n+1;i++)  */
	/* 	for(j=1;j<n+1;j++) */
	/* 	  M2->P[i][j] = M->P[i][j];  */

	/*       //print_int_matrix(M->P,n+1,n+1); */
	/*       //print_int_matrix(M2->P,n+1,n+1); */


	/*       /\* Choose random edge *\/ */
	/*       R1 = (int) R4[iter] / n + 1; */
	/*       R2 = R4[iter] % n + 1; */
	/*       if (R1 == R2) */
	/* 	reject = 1; */
	/*       iter++; */

	/*       //      printf("R1=%i\tR2=%i\n",R1,R2); */

	/*       if (!reject) */
	/* 	{ */
	/* 	  M2->P[R1][R2] = 1 - M2->P[R1][R2]; // Change edge */

	/* 	  if (M2->P[R1][R2]) // added edge (needs testing) */
	/* 	    { */
	/* 	      poset_stat = reduce_to_cover_relations (M2->P, n); */
	/* 	      if (poset_stat != 2) // Fatal */
	/* 		{ */
	/* 		  if (!M2->P[R1][R2]) // ie, suggested relation has been removed */
	/* 		    { */
	/* 		      c = 0; */
	/* 		      for(i=1;i<n+1;i++) */
	/* 			if(M2->P[R1][i] && M2->P[i][R2]) // Check if single intermediate state, ie R1 -> j -> R2 */
	/* 			  { */
	/* 			    j = i; */
	/* 			    c++; */
	/* 			  } */
	/* 		      if(c == 1) */
	/* 			{ */
	/* 			  M2->P[j][R2] = 0; // Remove j -> R2 */
	/* 			  M2->P[R1][R2] = 1; // Add R1 -> R2 */
	/* 			  if (reduce_to_cover_relations (M2->P, n) > 0) */
	/* 			    reject = 1; */
	/* 			} */
	/* 		      else reject = 1; */
	/* 		    } */
	/* 		  compatibility(D, N_u, M2); */
	/* 		  N_compatible = 0; */
	/* 		  for (k=0; k<N_u; k++) */
	/* 		    { */
	/* 		      N_compatible += D[k].is_compatible * D[k].count; */
	/* 		    } */
	/* 		  alpha_new = (double) N_compatible / (double) N; */
	/* 		  if (alpha - alpha_new > 0.2) */
	/* 		    { */
	/* 		      loglik_next[R1-1][R2-1] = -1/0.0; */
	/* 		      //printf("%i %i; d_alpha=%f\n",R1,R2,alpha-alpha_new); */
	/* 		      reject = 1; */
	/* 		    } */
	/* 		  //else printf("%i %i; D_alpha=%f\n",R1,R2,alpha-alpha_new); */
	/* 		} */
	/* 	      else reject = 1; */
	/* 	    } */
	/* 	} */

	iter = 0;

	/* Now chose edge and compute loglik */
	while (iter < n*n)
	{
		reject = 0;

		/* Copy poset*/
		for (i=1;i<n+1;i++)
			for(j=1;j<n+1;j++)
				M2->P[i][j] = M->P[i][j];

		/* Choose random edge */
		R1 = (int) R4[iter] / n + 1;
		R2 = R4[iter] % n + 1;

		/* Increases unexplained data?*/
		/* This is used as a heuristic to prevent very bad steps*/
		if (alpha_all[R1-1][R2-1] < alpha)
		{
			boltz = exp( (alpha_all[R1-1][R2-1] - alpha ) / 0.05 );
			R3 = (double) rand() / (double) RAND_MAX;
			if (R3 > boltz)
				reject = 1;
		}

		if(!reject)
		{

			M2->P[R1][R2] = 1 - M2->P[R1][R2]; // Change edge

			if (M2->P[R1][R2]) // added edge (needs testing)
			{
				poset_stat = reduce_to_cover_relations (M2->P, n);
				if (poset_stat != 2) // Fatal
				{
					if (!M2->P[R1][R2]) // ie, suggested relation has been removed
					{
						c = 0;
						for(i=1;i<n+1;i++)
							if(M2->P[R1][i] && M2->P[i][R2]) // Check if single intermediate state, ie R1 -> j -> R2
							{
								j = i;
								c++;
							}
						if(c == 1)
						{
							M2->P[j][R2] = 0; // Remove j -> R2
							M2->P[R1][R2] = 1; // Add R1 -> R2
							if (reduce_to_cover_relations (M2->P, n) > 0)
								reject = 1;
						}
						else reject = 1;
					}
				}
				else
				{
					reject = 1;
				}
			}

		}

		if (!reject)
		{
			printf("Testing:\n");
			//print_int_matrix(M2->P,n+1,n+1);
			printf("%i %i\n",R1,R2);

			make_model(M2, n, M2->P);
			*epsilon = 1- alpha_all[R1-1][R2-1]; // Set initial value for epsilon

			if (loglik_next[R1-1][R2-1] == 0.0) // Check if loglik has been computed before
			{
				guess_lambda(M2->P, n, D, N_u, lambda_new);
				//for(i=1;i<n+1;i++)
				//lambda_new[i] = lambda[i];
				loglik_next[R1-1][R2-1] = EM_EM(M2, D, N_u, lambda_new, epsilon);
			}
			loglik_new = loglik_next[R1-1][R2-1];

			printf("%f\t%f\n",loglik,loglik_new);

			if ( loglik_new > loglik) // Always accept increasing Loglik
			{
				make_model(M, n, M2->P);
				compatibility(D, N_u, M);
				for(i=1;i<n+1;i++)
					lambda[i] = lambda_new[i];
				free(R4);
				free(R5);
				for(i=0; i < M->n; i++)
					for(j=0; j < M->n; j++)
						loglik_next[i][j] = 0.0; // Re-initialize possible logliks
				return loglik_new;
			}
			else // Accept decreasing Loglik with Boltzmann weight
			{
				boltz = exp( (loglik_new - loglik ) / T );
				R3 = (double) rand() / (double) RAND_MAX;
				printf(" e^{-dL/T}=%f\tR=%f\n", boltz, R3);
				if ( R3 < boltz )
				{
					make_model(M, n, M2->P);
					compatibility(D, N_u, M2);
					for(i=1;i<n+1;i++)
						lambda[i] = lambda_new[i];
					free(R4);
					free(R5);
					for(i=0; i < M->n; i++)
						for(j=0; j < M->n; j++)
							loglik_next[i][j] = 0.0; // Re-initialize possible logliks
					return loglik_new;
				}
				else break;
			}
		}

		iter++;
	}
	printf("rejected\n");
	free(R4);
	free(R5);
	free(lambda_new);
	free(alpha_all);
	return loglik;
}

double** read_loglik_next(char* filestem, int n)
{
	int j, k;

	char suffix[15] = ".log";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	double** loglik_next = get_double_matrix(n,n);

	FILE *input;
	if ( (input = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "Warning:  Could not read %s\n", filename);
		for(k=0; k < n; k++)
			for(j=0; j < n; j++)
				loglik_next[k][j] = 0.0;
		return loglik_next;
	}

	/* Read patterns */
	float x;
	for (k=0; k<n; k++)
	{
		for (j=0; j<n; j++)
		{
			if (fscanf(input,"%g\t", &x) == 1)
			{
				if (x > 0)
				{
					fprintf(stderr, "ERROR reading data from %s!\n", filename);
					exit(1);
				}
				loglik_next[k][j] = (double)x;
			}
			else
			{
				fprintf(stderr, "Error reading data from %s!\n", filename);
				exit(1);
			}

		}
	}

	return loglik_next;
}

void write_loglik_next(char* filestem, double** loglik_next, const int n)
{
	int j,k;

	char suffix[15] = ".log";
	char *filename = (char *) calloc(strlen(filestem) + strlen(suffix) + 1, sizeof(char));
	strcat(filename, filestem);
	strcat(filename, suffix);

	FILE *output;
	if ( (output = fopen(filename, "w")) == NULL)
	{
		fprintf(stderr, "Error:  Could not write to file %s\n", filename);
		exit(1);
	}

	for (k=0; k<n; k++)
	{
		for (j=0; j<n; j++)
		{
			fprintf(output, "%lf\t", loglik_next[k][j]);
		}
		fprintf(output, "\n");
	}

	fclose(output);
}


double local_search(model* M, data* D, int N_u, double* lambda, double* epsilon, double loglik, double T, int N_iter, char* filestem, int l_flag, int w_flag)
{
	model M2;
	make_model(&M2, M->n, M->P);

	int iter = 0;

	double loglik_new = loglik;
	double** loglik_next = get_double_matrix(M->n,M->n); // Store loglik of possible steps
	int i,j;
	if(l_flag)
		loglik_next = read_loglik_next(filestem, M->n);
	else
	{
		for(i=0; i < M->n; i++)
			for(j=0; j < M->n; j++)
				loglik_next[i][j] = 0.0;
	}
	//print_double_matrix(loglik_next, M->n, M->n);

	printf("Step %i/%i\n===\n%f\t%f\n", iter,N_iter,T, loglik_new);
	//print_int_matrix(M->P,M->n+1,M->n+1);
	for(i=1;i<=M->n;i++)
		for(j=1;j<=M->n;j++)
			if(M->P[i][j])
				printf("%i %i\n",i,j);

	while (iter <= N_iter)
	{
		loglik_new = try_edge(M, &M2, D, N_u, lambda, epsilon, loglik_new, loglik_next, T);

		printf("===\nStep %i/%i\nT=%f\tL=%f\n", iter, N_iter, T, loglik_new);
		//print_int_matrix(M->P,M->n+1,M->n+1);
		for(i=1;i<=M->n;i++)
			for(j=1;j<=M->n;j++)
				if(M->P[i][j])
					printf("%i %i\n",i,j);
		printf("---\n");

		// Write output
		if (w_flag){
			write_poset(0, filestem, M->P, M->n, -1);
			write_loglik_next(filestem, loglik_next, M->n);
		}

		iter++;
		T *= 1 - 1 / (double)M->n;
	}

	printf("---\nLoglik=\n");
	print_double_matrix(loglik_next, M->n, M->n);
	for(i=0; i < M->n; i++)
		free(loglik_next[i]);
	free(loglik_next);

	return loglik_new;
}




void resample(data* D, double* p, int N_u)
{
	int k;
	unsigned int N = 0;
	unsigned int *n = get_uint_array(N_u);

	for (k=0; k<N_u; k++)
		N += D[k].count;

	gsl_ran_multinomial(RNG, N_u, N, p, n);

	for (k=0; k<N_u; k++)
		D[k].count = n[k];

	free(n);
}



int is_equal_int_matrix(int** A, int** B, int n)
{
	int i, j;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			if (A[i][j] != B[i][j])
				return 0;

	return 1;
}



void boolean_matrix_sum(int** A, int** B, int** C, int n)
{
	/*
      Boolean matrix sum  A + B = C
	 */

	int i, j;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			C[i][j] = A[i][j] + B[i][j];
			C[i][j] = C[i][j] ? 1 : 0;
		}

}



void int_matrix_sum(int** A, int** B, int** C, int n)
{
	/*
      Matrix sum  A + B = C
	 */

	int i, j;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			C[i][j] = A[i][j] + B[i][j];

}



void boolean_matrix_product(int** A, int** B, int** C, int n)
{
	/*
      Boolean matrix product  A * B = C
	 */

	int i, j, k;

	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
		{
			C[i][j] = 0;
			for (k=0; k<n; k++)
				C[i][j] += A[i][k] * B[k][j];
			C[i][j] = C[i][j] ? 1 : 0;
		}

}



void transitive_closure(int** A, int**T, int n)
{
	/*
      T is the transitive closure of the relation A
	 */

	int** R = get_int_matrix(n, n);
	int** S = get_int_matrix(n, n);

	int i, j, k = 0;

	// initialize T = A
	for (i=0; i<n; i++)
		for (j=0; j<n; j++)
			T[i][j] = A[i][j];

	while ((k == 0) || (! is_equal_int_matrix(S, T, n)))
	{
		// S = T
		for (i=0; i<n; i++)
			for (j=0; j<n; j++)
				S[i][j] = T[i][j];

		// T = A*S + A
		boolean_matrix_product(A, S, R, n);
		boolean_matrix_sum(R, A, T, n);

		k++;
	}

	for (i=0; i<n; i++)
	{
		free(R[i]);
		free(S[i]);
	}
	free(R);
	free(S);

}


int ML_path(model *M, double *lambda){
	int i,j,c;
	double lambda_next;
	int mut_next, index_next;

	i = 0;
	while(M->N_ch[i] > 0 && !(M->N_ch[i]==1 && M->ch_diff[i][0]==0))
	{
		lambda_next = 0.0;
		mut_next = 0;
		index_next = 0;
		for (c=0; c < M->N_ch[i]; c++) // loop
		{
			j = M->ch_diff[i][c];  // index of differing event
			if(lambda[j] > lambda_next && j > 0){
				mut_next = j;
				index_next = M->ch[i][c];
				lambda_next = lambda[j];
			}
		}
		i = index_next;
		printf("%i ", mut_next);
	}
	printf("\n");
	return 0;
}
