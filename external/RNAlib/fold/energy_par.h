#include "common.h"
/* 
   prototypes for energy_par.c
*/

#include "energy_const.h"

CPPEXTERN double lxc37;   /* parameter for logarithmic loop
			  energy extrapolation            */

CPPEXTERN int stack37[NBPAIRS+1][NBPAIRS+1];
CPPEXTERN int enthalpies[NBPAIRS+1][NBPAIRS+1]; /* stack enthalpies */
CPPEXTERN int entropies[NBPAIRS+1][NBPAIRS+1];  /* not used anymore */

CPPEXTERN int hairpin37[31];
CPPEXTERN int bulge37[31];
CPPEXTERN int internal_loop37[31];
CPPEXTERN int internal2_energy;
CPPEXTERN int old_mismatch_37[NBPAIRS+1][5][5];
CPPEXTERN int mismatchI37[NBPAIRS+1][5][5];  /* interior loop mismatches */
CPPEXTERN int mismatchH37[NBPAIRS+1][5][5];  /* same for hairpins */
CPPEXTERN int mismatchM37[NBPAIRS+1][5][5];  /* same for multiloops */
CPPEXTERN int mism_H[NBPAIRS+1][5][5];       /* mismatch enthalpies */

CPPEXTERN int dangle5_37[NBPAIRS+1][5];      /* 5' dangle exterior of pair */
CPPEXTERN int dangle3_37[NBPAIRS+1][5];      /* 3' dangle */
CPPEXTERN int dangle3_H[NBPAIRS+1][5];       /* corresponding enthalpies */
CPPEXTERN int dangle5_H[NBPAIRS+1][5];

CPPEXTERN int int11_37[NBPAIRS+1][NBPAIRS+1][5][5]; /* 1x1 interior loops */
CPPEXTERN int int11_H[NBPAIRS+1][NBPAIRS+1][5][5];

CPPEXTERN int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5]; /* 2x1 interior loops */
CPPEXTERN int int21_H[NBPAIRS+1][NBPAIRS+1][5][5][5];

CPPEXTERN int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5]; /* 2x2 interior loops */
CPPEXTERN int int22_H[NBPAIRS+1][NBPAIRS+1][5][5][5][5];

/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */
CPPEXTERN int ML_BASE37;
CPPEXTERN int ML_closing37;
CPPEXTERN int ML_intern37;

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
CPPEXTERN int         MAX_NINIO;                   /* maximum correction */
CPPEXTERN int F_ninio37[5];

/* penalty for helices terminated by AU (actually not GC) */
CPPEXTERN int TerminalAU;
/* penalty for forming bi-molecular duplex */
CPPEXTERN int DuplexInit;
/* stabilizing contribution due to special hairpins of size 4 (tetraloops) */
CPPEXTERN char Tetraloops[];  /* string containing the special tetraloops */
CPPEXTERN int  TETRA_ENERGY37[];  /* Bonus energy for special tetraloops */
CPPEXTERN int  TETRA_ENTH37;
CPPEXTERN char Triloops[];    /* string containing the special triloops */
CPPEXTERN int  Triloop_E37[]; /* Bonus energy for special Triloops */  

CPPEXTERN double Tmeasure;       /* temperature of param measurements */
