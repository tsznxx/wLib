#include "common.h"
CPPEXTERN float  alifold(char **strings, char *structure);
CPPEXTERN  void  free_alifold_arrays(void);
CPPEXTERN double cv_fact /* =1 */;
CPPEXTERN double nc_fact /* =1 */;
typedef struct {
   short i;        /* i,j in [0, n-1] */
   short j;
   float p;      /* probability */
   float ent;    /* pseudo entropy for p(i,j) = S_i + S_j - p_ij*ln(p_ij) */
   short bp[8];  /* frequencies of pair_types */
   char comp;    /* 1 iff pair is in mfe structure */
}  pair_info;

/* CPPEXTERN float aliLfold(char **strings, char *structure, int maxdist); */
CPPEXTERN float alipf_fold(char **sequences, char *structure, struct plist **pl);
/* CPPEXTERN float alipfW_fold(char **sequences, char *structure, struct plist **pl,int winSize, float cutoff, int pairsize); */
/* CPPEXTERN struct plist *get_plistW(struct plist *pl, int length, double cutoff, int start, FLT_OR_DBL **Tpr, int winSize); */
CPPEXTERN char *centroid_ali(int length, double *dist,struct plist *pl);
CPPEXTERN float **readribosum(char *name);
CPPEXTERN char *alipbacktrack(double *prob) ;
CPPEXTERN void  free_alipf_arrays(void);
CPPEXTERN float  energy_of_alistruct(char **sequences, const char *structure, int n_seq, float *CVenergy);
CPPEXTERN float circalifold(const char **strings, char *structure);
CPPEXTERN float alipf_circ_fold(char **sequences, char *structure, struct plist **pl);
