#include <string>
#include <stdlib.h>
#include <vector>
#if !defined(STRUCTURE_H)
#define STRUCTURE_H




#include "defines.h"

#ifdef EXTENDED_DOUBLE
	#include "extended_double.h" //inlcude code for extended double if needed
#endif//defined EXTENDED_DOUBLE

//////////////////////////////////////////////////////////////////////
class structure //this structure contains all the info for a structure
{
public:
int numofstructures;//number of structures 
int numofbases;//number of nucleotides in sequence
short int pair[maxforce][2],npair,nforbid,forbid[maxforce][2];//arrays to hold lists of forced pairs or pairs forbidden
short int *numseq,*hnumber;
int **basepr;
int *energy;//[maxstructures+1];
char **ctlabel;//[maxstructures+1][ctheaderlength];
short int ndbl, dbl[maxforce];
int inter[3],allocatedstructures;
short int nnopair,*nopair,nmod,mod[maxforce];
int nopairmax;//maximum number of nucleotides not allowed to pair
short int ngu,gu[maxgu];
char *nucs;
bool intermolecular,allocated,templated,stacking;
bool **tem;//tem stores template information as to whether a pair is allowed
bool SHAPEFileRead,distsread;
structure(int structures = maxstructures+1);
~structure();
void allocate(int size = maxbases);
void allocatestructure();
void deletestructure();
void checknopair();//make sure that the nopair array is large enough to add more items
void checknumberofstructures();//check to make sure there is room for one more structure
void allocatetem();//allocate space in **tem 
short int min_gu, min_g_or_u;//NMR-derived constraint variables
short int neighbors[maxforce][maxneighborlength],nneighbors;//also NMR-derived index this from zero in both dimensions
//regional NMR constraints:
short int nregion,rmin_gu[maxregions],rmin_g_or_u[maxregions];
short int rneighbors[maxregions][maxforce][maxneighborlength],rnneighbors[maxregions],start[maxregions],stop[maxregions];
//microarray type constraints:
short int nmicroarray,microstart[maxregions],microstop[maxregions],microunpair[maxregions];
bool *fcedbl;//pointer to a 2-D array used in Dynalign to track nucleotides that must be double-stranded
void sort();//sort the structures so that the lowest free energy structure is in position one
			//NOTE: This function sorts the the list of base pairs according to energy.  The energies (*energy) and basepairs (**basepr)
			//		arrays are correctly sorted, but ct labels (**ctlabel) is not sorted.
double CalculatePseudoEnergy(double data, std::string modifier, double, double);
double Gammadist(double data, double shape, double loc, double scale);
double Potential(double data, std::vector< std::vector<double> > params, double kT);
void ReadProbabilisticPotentialParams();//Read chemical modifier distributions from file
void ReadSHAPE(const char *filename, std::string modifier="SHAPE", bool calculate=true, bool nosum=false);//Read SHAPE reactivity data from a file
void ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold);//Read SHAPE reactivity data from a file
void ReadOffset(const char *SSOffset, const char *DSOffset);//Read Free Energy Offset Files.
void ReadExperimentalPairBonus(const char *filename, double const experimentalOffset = 0.0, double const experimentalScaling = 1.0 );
bool limitdistance;//toggle to indicate that there is a limit on the maximum distance between nucs in base pairs
int maxdistance;//maximum distance between nucs in base pairs
double *SHAPE;//double array to contain SHAPE data -- values less than -500 are ignored
double **EX;// double array that contains experimental bonuses/penalties
bool shaped;//keeps track of whether SHAPE data was loaded
bool experimentalPairBonusExists;//keeps track of whether experimental bonus data was loaded
bool ssoffset;//keeps track of whether a single stranded offset was read from disk
double SHAPEslope,SHAPEintercept;//values of slope and intercept for SHAPE data modification of pairing stability
//SINGLE STRANDED SHAPE ENERGY VARIABLES AND FUNCTIONS
double *SHAPEss; //short int array that contains SHAPE data for single-stranded segments
double SHAPEslope_ss, SHAPEintercept_ss; //values of the slope and intercept for SHAPE data modifying single stranded loop stability
short int **SHAPEss_region;  //2-d short int array containing energy values for hairpin loop combinations
//Parameters for distributions
std::vector< std::vector<double> > SHAPE_params;
std::vector< std::vector<double> > DMS_params;
std::vector< std::vector<double> > CMCT_params;
int SHAPEss_calc(int index_i, int index_j);  //Returns pseudoenergy term for a hairpin loop using single stranded SHAPE data
short int SHAPEss_give_value(int index);  //Returns the single stranded SHAPE pseudo energy for a given nucleotide

void writedotbracket(const char *filename);//write dot-bracket notation of ct file (no pseudoknots!!)

double **constant;//constant is used to hold an array of equilibrium constants.  In partition function calculations, 
					//the equilibrium constant is multiplied by constant[j][i] when the i-j pair is formed. 
					//NOTE: The use of constant is NOT orthogonal to using chemical modification data.  They cannot
					//both be used at once.
void allocateconstant();//Function to allocate memory for constant array.

/*structure is set up to hold many possible structures of the same sequence
	numofbases = number of bases in sequence
	numofstructures = number of alternative structures of the sequence
				that is held by structure
	numseq[i] = a numeric that stands for the base in the ith position
				of the sequence,
			A = 1
			C = 2
			G = 3
			U = 4
	basepr[i][j] = base to which the jth base is paired in the ith structure
	force[i] = any information about the way the ith pair is treated when
				folding; eg: forced single, etc
	energy[i] = 10 * the Gibb's free energy of the ith structure, this
				is done so that integer math can be performed
				and the nearest tenth of a kcal/mol can be
				followed
	ctlabel = a string of information for each of the structures
	fce = an array that can keep track of how each i,j pair is being
				forced
   hnumber array stores the historical numbering of a sequence
   nucs is a character array to store the sequence information --
   	this will allow the program to keep T and U from getting confused

	stacking = a bool to keep track of whether stacking is being tracked in the ct.
		If stacking is to be tracked, stacking must be manually set to true before allocate
		is called.  The stacks are stored in basepr, by doubling the size of basepr, in basepr[index][i+N],
		where N is the number of nucleotides.  This lets the stacking info get sorted along with the
		pairing info in ::sort.
		When stacking is true, basepr[i][j]!=implies a stack for nucleotide j-numofbases.
		For j-numofbases unpaired and j-numofbases is not indicated in another stack, this nucleotide (j-numofbases) 
		is either half a terminal mismatch or a dangling end that is stacked on the indicated nucleotide.
		For j-numofbases paired, this is a coaxial stack in which j-numofbases is stacked onto the 
		nucleotide as indicated.  For flush stacks, both nucleotides in terminal basepairs indicate each other as
		stacks.  For stacks with an intervening mismatch, one helix will point to an unpaired nuc that mediates the stack
		and then that unpaired nucleotide will indicate the next helix in the coaxial stack.

	
*/



};

short int tonumi(char *base); //converts base to a numeric


char *tobase (int i);//convert a numeric value for a base to the familiar
								//character
void sortstructures (structure *ct);//this routine resorts the structures
					//predicted by dynamic according to energies calculated by efn2

void tonum(char *base,structure *ct,int count); //converts base to a numeric

void openct(structure *ct,const char *ctfile);//reads ct file

int openseq (structure *ct,const char *seqfile);//inputs a sequence from file
															// seqfile

//outputs a ct file, provide pointer to structure, pointer to cstring with filename
//	if append is set to true, the ct information is appended to an existing file, if the file exists.
void ctout (structure *ct,const char *ctoutfile, bool append=false);


int ecompare(const void *i, const void *j);

void swap(int *a,int *b);//Swap two variables
void swap(short int *a,short int *b);//swap two variables
void swap(float *a,float *b);//swap two variables
void swap(double *a,double *b);//swap two variables
#endif
