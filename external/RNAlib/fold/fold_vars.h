#include "common.h"
#ifndef __FOLD_VARS__
#define __FOLD_VARS__

/* to use floats instead of doubles in pf_fold() comment next line */
#define LARGE_PF
#ifdef  LARGE_PF
#define FLT_OR_DBL double
#else
#define FLT_OR_DBL float
#endif

CPPEXTERN int  noGU;           /* GU not allowed at all */
CPPEXTERN int  no_closingGU;   /* GU allowed only inside stacks */
CPPEXTERN int  tetra_loop;     /* Fold with specially stable 4-loops */
CPPEXTERN int  energy_set;     /* 0 = BP; 1=any mit GC; 2=any mit AU-parameter */
CPPEXTERN int  dangles;	    /* use dangling end energies (not in part_func!) */
/*@null@*/

CPPEXTERN int oldAliEn;        /* use old alifold energies (with gaps) */
CPPEXTERN int ribo;            /* use ribosum matrices */
CPPEXTERN char *RibosumFile;   /* warning this variable will vanish in the future
			       ribosums will be compiled in instead */
CPPEXTERN char *nonstandards;  /* contains allowed non standard bases */
CPPEXTERN double temperature;   /* rescale parameters to this temperature */
CPPEXTERN int  james_rule;     /* interior loops of size 2 get energy 0.8Kcal and
			       no mismatches, default 1 */
CPPEXTERN int  logML;          /* use logarithmic multiloop energy function */
CPPEXTERN int  cut_point;      /* first position of 2nd strand for co-folding */

typedef struct bond {               /* base pair */
   int i;
   int j;
} bondT;
CPPEXTERN bondT  *base_pair; /* list of base pairs */

CPPEXTERN FLT_OR_DBL *pr;          /* base pairing prob. matrix */
CPPEXTERN int   *iindx;            /* pr[i,j] -> pr[iindx[i]-j] */
CPPEXTERN double pf_scale;         /* scaling factor to avoid float overflows*/
CPPEXTERN int    fold_constrained; /* fold with constraints */
CPPEXTERN int    do_backtrack;     /* calculate pair prob matrix in part_func() */
CPPEXTERN int    noLonelyPairs;    /* avoid helices of length 1 */
CPPEXTERN char backtrack_type;     /* usually 'F'; 'C' require (1,N) to be bonded;
				   'M' seq is part of a multi loop */
char * option_string(void);
#endif
