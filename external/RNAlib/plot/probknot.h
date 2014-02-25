

#include "structure.h"
#include "pfunction.h"


//Build the ProbKnot structure.

//return an int that indicates errors.  0 = no error.
//This requires: v, a pointer to pfunction class, which is the v array of a previous partition function calculation.
//w5, a pointer to PFPRECISION, the w5 array of a previous partition function calculation.
//ct, a pointer to structure, which is filled with the sequence data and will receive the structure data.
//data, a pointer to pfdatatable, which has the thermodynamic parameters.
//lfce, a pointer a bool array, as filled by a previous partition function calculation.
//mod, a pointer a bool array, as filled by a previous partition function calculation.
//scaling, the scaling per nucleotide used by the partition function calculation.
//fce, a pointer to forceclass, as used by the previous partition function calculation.
//iteration, an int that indicates the number of assembly iteration to be performed, the defaulty and recommended value are 1.
//MinHelixLength, and int that indicates the shortest helix length allowed.  This defaulst to 1, the recommended value.
int ProbKnotAssemble(pfunctionclass *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce, int iterations =1, int MinHelixLength=1 );


//Remove short helices from a structure, allowing stacks across a single bulge
//Requires a pointer to ct that has the structure.
//Also requires the minimum helix length, MinHelixLength, and the number of the structure from which to remove pairs, StructureNumber, which defaults to 1.
void RemoveShortHelices(structure *ct, int MinHelixLength, int StructureNumber =1 );

