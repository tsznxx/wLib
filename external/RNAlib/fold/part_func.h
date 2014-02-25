#include "common.h"
#ifndef __PART_FUNC__
#define __PART_FUNC__

#define FLT_OR_DBL double
/* functions from part_func.c */
CPPEXTERN float  pf_fold(char *sequence, char *structure);
/* calculate partition function and base pair probabilities */
CPPEXTERN void   init_pf_fold(int length);    /* allocate space for pf_fold() */
CPPEXTERN void   free_pf_arrays(void);        /* free arrays from pf_fold() */
CPPEXTERN void   update_pf_params(int length); /*recalculate energy parameters */
CPPEXTERN char   bppm_symbol(const float *x);  /* string representation of structure */
CPPEXTERN double mean_bp_dist(int length); /* mean pair distance of ensemble */
CPPEXTERN char  *centroid(int length, double *dist);     /* mean pair distance of ensemble */
CPPEXTERN int get_pf_arrays(short **S_p, short **S1_p, char **ptype_p, FLT_OR_DBL **qb_p, FLT_OR_DBL **qm_p, FLT_OR_DBL **q1k_p, FLT_OR_DBL **qln_p);

CPPEXTERN	char	*pbacktrack(char *sequence);
CPPEXTERN	float	pf_circ_fold(char *sequence, char *structure);
CPPEXTERN	char	*pbacktrack_circ(char *seq);

#endif
