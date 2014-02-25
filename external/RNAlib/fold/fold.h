#include "common.h"
/* function from fold.c */
CPPEXTERN float  fold(const char *sequence, char *structure);
/* calculate mfe-structure of sequence */
CPPEXTERN float  energy_of_struct(const char *string, const char *structure);
/* calculate energy of string on structure */
CPPEXTERN void   free_arrays(void);           /* free arrays for mfe folding */
CPPEXTERN void   initialize_fold(int length); /* allocate arrays for folding */
CPPEXTERN void   update_fold_params(void);    /* recalculate parameters */
CPPEXTERN char  *backtrack_fold_from_pair(char *sequence, int i, int j);
CPPEXTERN int loop_energy(short * ptable, short *s, short *s1, int i);
CPPEXTERN void		export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);

/* some circfold related functions...	*/
CPPEXTERN	float	circfold(const char *string, char *structure);
CPPEXTERN	float	energy_of_circ_struct(const char *string, const char *structure);
CPPEXTERN	void	export_circfold_arrays(int *Fc_p, int *FcH_p, int *FcI_p, int *FcM_p, int **fM2_p, int **f5_p, int **c_p, int **fML_p, int **fM1_p, int **indx_p, char **ptype_p);
