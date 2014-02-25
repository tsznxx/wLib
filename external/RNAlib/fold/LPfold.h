#include "common.h"
#ifndef LPFOLD_H
#define LPFOLD_H

CPPEXTERN void  update_pf_paramsLP(int length);
CPPEXTERN struct plist *pfl_fold(char *sequence, int winSize, int pairSize, float cutoffb, double **pU, struct plist **dpp2, FILE *pUfp, FILE *spup);
CPPEXTERN void putoutpU_prob(double **pU,int length, int ulength, FILE *fp, int energies);
#endif
