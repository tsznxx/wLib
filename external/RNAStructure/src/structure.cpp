


//#if !defined(STRUCTURE_CPP)
//#define STRUCTURE_CPP

//#include <stdlib.h>
#include "structure.h"
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <math.h>
using namespace std;


structure::structure(int structures)
{
	allocatedstructures = structures;
	int i;

	allocatestructure();


	allocated = false;

	nnopair=0;
	nopairmax=maxforce-1;
	nopair= new short int [maxforce];//dynamically allocate nopair so that the list size can be expanded if need be
	nforbid=0;
	npair=0;
	ndbl=0;
	intermolecular = false;
	ngu = 0;
	templated = false;
	nmod = 0;
	min_gu = 0;
	min_g_or_u = 0;
	nneighbors = 0;
	nregion = 0;
	nmicroarray=0;
	numofstructures=0;

	//initialize values for SHAPE slope and intercept as zero for both double and single stranded restraints.
	SHAPEslope_ss = 0;
	SHAPEintercept_ss = 0;

	SHAPEslope = 0;
	SHAPEintercept = 0;


	for (i=0;i<maxregions;i++) rnneighbors[i]=0;
	//for (i = 0; i < 100;i++) {
	//	for ( j = 0; j < 25; j++) {
	//		neighbors[i][j]=0;
	//	}
	//}

	stacking = false;
	limitdistance=false;//toogle to true to limit base pairing distance
	maxdistance=600;//default maximum distance between paired nucs
	shaped = false;//by default, a structure does not have SHAPE data associated with it
	SHAPEss_region = NULL;//by default, the SHAPEss_region does not need to be allocated
	ssoffset = false;//by default, a structure does not have a single stranded offset read from disk
	experimentalPairBonusExists = false;//by default, no pairwise bonuses provided
	constant = NULL;//by default, indicate that the array constant is not allocated
	SHAPEFileRead = false;
	distsread = false;

}

void structure::allocatestructure() {
	int i;
	energy = new int [allocatedstructures];
	ctlabel = new char *[allocatedstructures];
	for (i=0;i<allocatedstructures;i++) ctlabel[i]=new char [ctheaderlength];
	for (i=1;i<allocatedstructures;i++) {
		energy[i]=0;
	}

}


//make sure the nopair array (listing nucleotides that must be single-stranded) is long enough
	//so that more nucleotides can be added to the list.
void structure::checknopair() {
	short int *temp;
	int i;

	if (nnopair >= nopairmax) {
	nopairmax= nopairmax*2;
	temp = new short int [nopairmax];

	for (i=1;i<=nnopair;i++) {
		temp[i]=nopair[i];
	}

	delete [] nopair;
	nopair=temp;
	}

}

void structure::deletestructure() {
	int i;

	if (allocated) {
		for (i=0;i<allocatedstructures;i++) {
    		delete[] basepr[i];
   		}

		delete[] basepr;
	}

	delete[] energy;
	for (i=0;i<allocatedstructures;i++) delete[] ctlabel[i];
	delete[] ctlabel;
}

structure::~structure()
{
	int i;



	deletestructure();
	delete[] nopair;
	if (allocated) {
		delete[] numseq;

		delete[] hnumber;
		delete[] nucs;
	}
	if (templated) {
		for (i=0;i<=numofbases;i++) {
    		delete[] tem[i];
   		}

   		delete[] tem;
	}
	if (shaped) {
		delete[] SHAPE;
		//delete ss shape array
		delete[] SHAPEss;
		if (SHAPEss_region!=NULL) {
			for (int i = 1; i <= numofbases; i++) {
				delete[] SHAPEss_region[i];
			}
			delete[] SHAPEss_region;
		}
	}
	if ( experimentalPairBonusExists ) {
		delete[] EX;
	}
	if (constant!=NULL) {
		//delete the equilibrium constants
		for (i=0;i<=numofbases;i++) {
    		delete[] constant[i];
   		}

   		delete[] constant;

	}
}



void structure::allocate(int size)

{
	int i;

   numseq = new short int [2*size+1];
   hnumber = new short int [size+1];
   nucs = new char [size+2];
   basepr = new int *[allocatedstructures];
   for (i=0;i<allocatedstructures;i++) {
    	if (stacking) basepr[i] = new int [2*size+1];
		else basepr[i] = new int [size+1];
   }
   allocated = true;



}


//this allocates space in an array that is used for folding with phylogenetic data
//	tem == template for allowed and disallowed pairs
void structure::allocatetem()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   tem = new bool *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	tem[i] = new bool [i+1];
   }
   templated = true;

   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i;j<=numofbases;j++) {
    		tem[j][i] = true;
		}
   }

}

//this allocates space in an array that is used for applying an equilibrium constant for missing specific pairs
void structure::allocateconstant()

{
	int i,j;
	//Size = size;//save the size of the array so that the destructor can
   				//deallocate the space

   constant = new double *[numofbases+1];
   for (i=0;i<=numofbases;i++) {
    	constant[i] = new double [i+1];
   }


   //initialize all positions to true:
	for (i=0;i<=numofbases;i++) {
		for (j=i;j<=numofbases;j++) {
    		constant[j][i] = 1.0;
		}
   }

}

void structure::checknumberofstructures() {
	structure *ct2;
	int i,j;

	if (numofstructures>=allocatedstructures-1) {
		//there will be no room for the next structure, double the allocated space for structures
		ct2=new structure (allocatedstructures);
		ct2->stacking = stacking;
		ct2->allocate(numofbases);

		for (i=1;i<allocatedstructures-1;i++) {
			strcpy(ct2->ctlabel[i],ctlabel[i]);
			ct2->energy[i]=energy[i];

			for (j=1;j<=numofbases;j++) {
				ct2->basepr[i][j] = basepr[i][j];
				if (stacking) ct2->basepr[i][j+numofbases]=basepr[i][j+numofbases];
			}

		}


		deletestructure();

		allocatedstructures=max(2*allocatedstructures,numofstructures+1);

		allocatestructure();
		basepr = new int *[allocatedstructures];
		for (i=0;i<allocatedstructures;i++) {
    		if (stacking) basepr[i] = new int [2*numofbases+1];
			else basepr[i] = new int [numofbases+1];
		}


		for (i=1;i<=numofstructures-1;i++) {
			strcpy(ctlabel[i],ct2->ctlabel[i]);
			energy[i]=ct2->energy[i];

			for (j=1;j<=numofbases;j++) {
				basepr[i][j] = ct2->basepr[i][j];
				if (stacking) basepr[i][j+numofbases] =ct2->basepr[i][j+numofbases];
			}


		}

		delete ct2;
	}

}

//sort the structures from lowest to highest free energy
void structure::sort() {
	int i;

	for (i=1;i<=numofstructures;i++) basepr[i][0]=(short int) energy[i];

	//take advantage of qsort--
	//store the energies in basepr[i][0] so that they are associated with that array

	basepr++;

	qsort((void*) ((basepr)),(size_t) numofstructures,/*(numofbases+1)**/(sizeof(short int *)),ecompare);

	basepr--;
	for (i=1;i<=numofstructures;i++) energy[i]=(int)basepr[i][0];

}
double structure::Gammadist(double data, double shape, double loc, double scale){
	return (1/scale)*pow((data - loc)*(1/scale), (shape - 1))*exp(-(1/scale)*(data - loc))/tgamma(shape);
}


double structure::Potential(double data, std::vector< std::vector<double> > params, double kT){
	// params[0] is for paired, params[0] for unpaired...params[][j], j=0,1,2 for shape, loc scale of component 1
	// j=3,4,5 for shape, loc, scale of component 2 and j=6,7 for weights of components 1 and 2 respectively.
	double pairedprob = params[0][6]*Gammadist(data, params[0][0], params[0][1], params[0][2]) + 
	                   params[0][7]*Gammadist(data, params[0][3], params[0][4], params[0][5]); 
	double unpairedprob = params[1][6]*Gammadist(data, params[1][0], params[1][1], params[1][2]) + 
	                     params[1][7]*Gammadist(data, params[1][3], params[1][4], params[1][5]);
	return -kT*log(pairedprob/unpairedprob);
}

void structure::ReadProbabilisticPotentialParams() {
	string filedir;
	char *dir=getenv("DATAPATH");
	
	//Set filedir to DATAPATH, of DATAPATH exists as an environment variable:
	if (dir!=NULL) {
		filedir=dir;
	}
	//If DATAPATH does not exist, assume files are in pwd:
	else {
		filedir="";
	}
	filedir += "dists/";
	string line;
	int start, end, i, j;
	int nparams = 8;
	// Initialize parameters
	std::vector<double> shapecol1;
	std::vector<double> shapecol2;
	for(i=0; i<nparams; i++){
		shapecol1.push_back(0.0);
		shapecol2.push_back(0.0);
	}
	SHAPE_params.push_back(shapecol1);
	SHAPE_params.push_back(shapecol2);

	std::vector<double> dmscol1;
	std::vector<double> dmscol2;
	for(i=0; i<nparams; i++){
		dmscol1.push_back(0.0);
		dmscol2.push_back(0.0);
	}
	DMS_params.push_back(dmscol1);
	DMS_params.push_back(dmscol2);

	std::vector<double> cmctcol1;
	std::vector<double> cmctcol2;
	for(i=0; i<nparams; i++){
		cmctcol1.push_back(0.0);
		cmctcol2.push_back(0.0);
	}
	CMCT_params.push_back(cmctcol1);
	CMCT_params.push_back(cmctcol2);

	// Start reading distribution parameters....note that this code is ugly, can be collapsed into a single parameter array
	// SHAPE
	string tmp = filedir + "SHAPEdist.txt";
	char *filename = (char*)tmp.c_str();
	ifstream SHAPEfile;
	SHAPEfile.open(filename);
	if (SHAPEfile.is_open()){
		getline(SHAPEfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(SHAPEfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				SHAPE_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		SHAPEfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}

	// DMS
	tmp = filedir + "DMSdist.txt";
	filename = (char*)tmp.c_str();
	ifstream DMSfile;
	DMSfile.open(filename);
	if (DMSfile.is_open()){
		getline(DMSfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(DMSfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				DMS_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		DMSfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}

	// CMCT
	tmp = filedir + "CMCTdist.txt";
	filename = (char*)tmp.c_str();
	ifstream CMCTfile;
	CMCTfile.open(filename);
	if (CMCTfile.is_open()){
		getline(CMCTfile, line);
		for(i=0;i<2;i++) {
			start = 0;
			end = 0;
			getline(CMCTfile, line);
			for(j=0;j<nparams;j++) {
				end = line.find(" ", start);
				CMCT_params[i][j] = atof(line.substr(start, end).c_str());
				start = end + 1;
			}
		}
		CMCTfile.close();
	}
	else {
		cout << "Cannot open file " + tmp + "\n";
	}
}

// This function calculates the pseudoenergy for a given reactivity data. It changes the calculation
// depending on the modifier specified, giving either the log-likelihood-ratio of the unpaired/paired
// probabilities given a reactivity distribution per modifier, or the "classic" Deigan et al. bonus
// term when no modifier or an unrecognized modifier is provided. 
double structure::CalculatePseudoEnergy(double data, std::string modifier, double slope, double intercept) {
	std::vector< std::vector<double> > params;
	if( data <= -500)
		return 0;
	if(modifier == "SHAPE_AC" || modifier == "SHAPE_GU") {
		// This is only applied if SHAPE_AC or SHAPE_GU is specified
		// For now, I'm using the "default" calculations for SHAPE
		// pseudoenergies when the modifier is "SHAPE".
		params = SHAPE_params; 
	}
	else { 
		if( modifier == "DMS" ) {
			params = DMS_params; 
		}
		else{ 
			if( modifier == "CMCT") {
				params = CMCT_params; 
			}
			else {
				if (modifier == "diffSHAPE") {
					if (data>0){
						return data * slope;
					}
					else{
						return 0;
					}
				}
				else{
					if (data > 0) {
						return log( data + 1.0 )*slope + intercept;
					}
					else {
						return intercept;
						}
					}
				}
			}

		}
	if( data < 0 || (slope == 0 && intercept == 0))
		return 0;
	double val2 = log(data+1.0)*slope+intercept;
	double kT = 5.904976983149999;
	double val = Potential(data, params, kT);

#ifndef WIN32
	//This is only included if this is not compilig on Windows, where isnan is not available.
	if (isnan(val)!=0) {
		return 0;

	}
#endif
	if( (slope == 0 && intercept == 0) ){
		return 0;
	}
	else{
		return val;
	}

}


//This function reads a SHAPE reactivity datafile and saves the data for a linear penalty.
//calculate (default true) indicate whether these data are being read for folding.  (false means
	//the raw values need to be stored.)
void structure::ReadSHAPE(const char *filename, std::string modifier, bool calculate, bool nosum) {
	ifstream in;
	int position;
	double data;

	//Only read the dists if they need to be read.  These are not used for  SHAPE or diffSHAPE.
	if( !distsread && !(modifier=="SHAPE"||modifier=="diffSHAPE") ){

		ReadProbabilisticPotentialParams();
		distsread=true;
	}

	if( !SHAPEFileRead ){
		SHAPE = new double [2*numofbases+1];

		

		//initializes an array to hold SHAPE data for single stranded segments
		SHAPEss = new double [2 * numofbases + 1];

		SHAPEFileRead = true;

		for (position=0; position <= 2*numofbases; position++) {
			SHAPE[position] = 0;
			SHAPEss[position] = 0;
		}
	}

	double *SHAPEnew = new double [2*numofbases+1];

	int *counts = new int [numofbases+1];

	//initializes a temporal array to hold reactivity data for single stranded segments
	double *SHAPEssnew = new double [2 * numofbases + 1];

	// useful for bootstrapping -- keep count of how many times the data for a given position is specified.
	int *num_data_points =  new int [ numofbases + 1];
	for (position=0; position <= numofbases; position++) num_data_points[ position ] = 0;

	in.open(filename);
	shaped = true;

	for (position=0; position <= 2*numofbases; position++) {
		SHAPEnew[position] = 0.0;
		SHAPEssnew[position] = 0.0;
	}
	for (position=0; position <= numofbases; position++) {
		counts[position] = 0;
	}

	in >> position;
	in >> data;

	while (!in.eof()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity

		if (position<=numofbases) {
			if (calculate) { 
				SHAPEnew[position] += CalculatePseudoEnergy(data, modifier, SHAPEslope, SHAPEintercept);
				SHAPEssnew[position] += CalculatePseudoEnergy(data, modifier, SHAPEslope_ss, SHAPEintercept_ss);
				if(SHAPEnew[position] != 0)
					counts[position] += 1;
			}

			else {
				SHAPE[position] = data;
				SHAPEss[position] = data;
				num_data_points[ position ] = 1;
			}

		}


		in >> position;
		in >> data;
	}
	in.close();
	if (calculate) {
		for (position=1;position<=numofbases;position++) {
			if(counts[position] >= 1){
				SHAPE[position]+=SHAPEnew[position]/counts[position];
				SHAPEss[position]+=SHAPEssnew[position]/counts[position];
			}
		}

		//Add the new SHAPE data to existing -- This change was made so that differential SHAPE can be added to existing SHAPE data.
		for (position=1;position<=numofbases;position++) {
			SHAPE[position+numofbases]=SHAPE[position];
			SHAPEss[position+numofbases]=SHAPEss[position];
		}
		
	}

	
	
	/*
	if (calculate) {
		for (position=1;position<=numofbases;position++) {
			if (SHAPE[position]>0) {
				SHAPE[position]= (log(SHAPE[position]+1.0)*SHAPEslope+SHAPEintercept);
				SHAPE[position] *= num_data_points[ position ]; //for bootstrapping.
			}
			else if (SHAPE[position]>-500) {
				SHAPE[position]= SHAPEintercept;

			}
			else {
				SHAPE[position]=0.0;

			}
			SHAPE[position+numofbases]=SHAPE[position];

			//repeat process for single stranded SHAPE array
			//can be incorporated into the above code since the initial arrays are identical

			if (SHAPEss[position] > 0) {
				SHAPEss[position] = log(SHAPEss[position] + 1.0) * SHAPEslope_ss + SHAPEintercept_ss;
				SHAPE[position] *= num_data_points[ position ]; //for bootstrapping.
			}
			else if (SHAPEss[position] > -500) {
				SHAPEss[position] = SHAPEintercept_ss;

			}
			else {
				SHAPEss[position] = 0.0;

			}
			//copy results for duplicate
			SHAPEss[position + numofbases] = SHAPEss[position];


		}
	}*/
	//initializing triangular 2-d array that stores ss SHAPE energies for loops. 1st index is ending location, 2nd index is starting location
	SHAPEss_region = new short int *[numofbases + 1];
	for (int i = 1; i <= numofbases; i++) SHAPEss_region[i] = new short int [i];

	//fills 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters
	for (int j = 2; j <= numofbases; j++) {
		SHAPEss_region[j][j - 1] = (short int)(SHAPEss[j] + SHAPEss[j-1]); //sets energy for "zero sized loop".  Acts as starting value to add onto below
		for (int i = j - 2; i >= 1; i--) {
			SHAPEss_region[j][i] = SHAPEss_region[j][i + 1] + (short int)(SHAPEss[i]); //adds energy for additional loop element
		}
	}

	delete[] num_data_points;
	delete[] SHAPEnew;
	delete[] SHAPEssnew;
	delete[] counts;


}

//Read offset files
//The SSOffset provides free energies to add to nucleotides if they are single-stranded
//The DSOffset provides free energies to add to nucleotides if they are double-stranded
//This function must be called ofter reading SHAPE files, if read, because the same infrastructure is used.

//Either filename can be NULL, in which case that offset is not recorded.
void structure::ReadOffset(const char *SSOffset, const char *DSOffset) {
	int i,position;
	double data;
	ifstream in,in2;

	if (!shaped) {

		//initialize the arrays because SHAPE wasn't previously read
		SHAPE = new double [2*numofbases+1];
		SHAPEss = new double [2 * numofbases + 1];


		for (i=1;i<2*numofbases;i++) {
			SHAPE[i]=0;
			SHAPEss[i]=0;

		}

		//also initialize the ssSHAPE lookup table for regions
		SHAPEss_region = new short int *[numofbases + 1];
		for (int i = 1; i <= numofbases; i++) SHAPEss_region[i] = new short int [i];

		shaped = true;
	}


	if (SSOffset!=NULL) {

			ssoffset = true;//we are reading an offset

		in.open(SSOffset);

		in >> position;
		in >> data;

		while (!in.eof()) {
			//read and parse all data
				//required format is rows with sequence position followed by reactivity
				//multiply the specified free energy change by converstionfactor, the factor by which energies in kcal/mol are multipled

			if (position<=numofbases) {
				SHAPEss[position]+=(data*conversionfactor);
				SHAPEss[position+numofbases]+=(data*conversionfactor);
			}

			in >> position;
			in >> data;
		}
		in.close();
	}

	if (DSOffset!=NULL) {
		in2.open(DSOffset);

		in2 >> position;
		in2 >> data;

		while (!in2.eof()) {
			//read and parse all data
				//required format is rows with sequence position followed by reactivity
				//multiply the specified free energy change by converstionfactor, the factor by which energies in kcal/mol are multipled


			if (position<=numofbases) {
				SHAPE[position]+=(data*conversionfactor);
				SHAPE[position+numofbases]+=(data*conversionfactor);

			}

			in2 >> position;
			in2 >> data;
		}
		in2.close();
	}

	//fills 2-d array with pseudo energy terms for loops from i-j using ss SHAPE parameters or offsets
	//If SHAPE was previously read, this is a redo of the action
	for (int j = 2; j <= numofbases; j++) {
		SHAPEss_region[j][j - 1] = (short int)(SHAPEss[j] + SHAPEss[j-1]); //sets energy for "zero sized loop".  Acts as starting value to add onto below
		for (int i = j - 2; i >= 1; i--) {
			SHAPEss_region[j][i] = SHAPEss_region[j][i + 1] + (short int)(SHAPEss[i]); //adds energy for additional loop element
		}
	}



}

//this function returns a psuedo energy term for a single stranded nucleotide based on SHAPE data
short int structure::SHAPEss_give_value(int index) {
	if (shaped) {
		if (index > numofbases) return (short int)SHAPEss[index - numofbases];
		else return (short int)SHAPEss[index];
	}
	else return 0;
}

//this is the function that will return a pseudo energy term for hairpin loops based off of SHAPE data
int structure::SHAPEss_calc(int index_i, int index_j) {
	if (shaped) {
		//accounts for the SHAPEss_region array being only NxN, not 2Nx2N
		if (index_i > numofbases) index_i -= numofbases;
		if (index_j > numofbases) index_j -= numofbases;
		if (index_i > index_j) {
			int temp_index = index_i;
			index_i = index_j;
			index_j = temp_index;
		}
		return SHAPEss_region[index_j][index_i];
	} else return 0;  //if no shaped data is being used, return zero
}


//This function reads a SHAPE reactivity datafile and parse the data into single-stranded amd chemical modification constraints.
void structure::ReadSHAPE(const char *filename, float SingleStrandThreshold, float ModificationThreshold) {
	ifstream in;
	int position;
	float data;

	in.open(filename);
	in >> position;
	in >> data;

	while (!in.eof()) {
		//read and parse all data
			//required format is rows with sequence position followed by reactivity

		if (position<=numofbases) {
			if (data>=SingleStrandThreshold) {
				nnopair++;
				nopair[nnopair]=position;
			}
			else if (data>=ModificationThreshold) {
				nmod++;
				mod[nmod]=position;
			}
		}
		in >> position;
		in >> data;
	}
	in.close();

}

// This function reads an experimental pair bonus file, similar to SHAPE, but just straightforward
// application as kcal bonuses.  As with SHAPE, bonus is applied at 0x, 1x, and 2x for
//  single stranded, edge base pairs, and internal base pairs.
void structure::ReadExperimentalPairBonus(const char *filename, double const experimentalOffset, double const experimentalScaling ) {

	//int position;
	//float data;

	//initializing 2-d array that stores bonuses for base pairs.
	// actually this is already done in RNA.cpp. Thanks Pablo!
	//delete[] EX;

    	ifstream in(filename);


	EX = new double *[ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++) EX[i] = new double [ 2*numofbases+1 ];
	for (int i = 0; i < 2*numofbases+1; i++)
	  for (int j = 0; j < 2*numofbases+1; j++)
	    EX[i][j] = 0.0;

	for (int i = 1; i <= 2*numofbases; i++)
	  for (int j = 1; j <= 2*numofbases; j++)
	    EX[i][j] = experimentalOffset * PFPRECISION( conversionfactor );
	experimentalPairBonusExists = true;

	int i( 1 ), j( 1 ), count( 0 );
	double val;

	if ( filename != "" ) {


	  while (!in.eof() && j <= numofbases) {

	    //read and parse all data
	    //required format is bonuses in square matrix
	    in >> val;

	    EX[i           ][j           ] += val * PFPRECISION( conversionfactor ) * experimentalScaling;
	    EX[i+numofbases][j           ] = EX[i][j];
	    EX[i           ][j+numofbases] = EX[i][j];
	    EX[i+numofbases][j+numofbases] = EX[i][j];
	    count++;

	    i++;
	    if ( i > numofbases ){
	      i = 1;
	      j++;
	    }

	  }
	  in.close();

	  if ( count != numofbases * numofbases ){
	    cerr << "Found too few values in experimental bonus file, found "<< count <<" but expected " << numofbases*numofbases << endl;
	    exit(1);
	  }
	}

}


//write a dot-bracket file
//Note:  This function assumes that there are no pseudoknots, which would make the output un-parsable

void structure::writedotbracket(const char *filename) {
	ofstream out;
	int i,j;


	out.open(filename);



	for (i=1;i<=numofstructures;i++) {



		out << "> " << ctlabel[i] << "\n";
		for (j=1;j<=numofbases;j++) {
			out << nucs[j];
		}
		out << "\n";
		for (j=1;j<=numofbases;j++) {
			if (basepr[i][j]>j) out << "(";
			else if (basepr[i][j]==0) out << ".";
			else out << ")";
		}
		out << "\n";


	}



	out.close();


}


//comparison function used by qsort in structure::sort
int ecompare(const void *i, const void *j) {

	return (**((short int**)i)-**((short int**)j));

}


//read a ct file with sequence and structural information
#define linelength 20



void openct(structure *ct,const char *ctfile) {
	int count,i,j;
	char base[2],header[ctheaderlength],temp[1000];
	istream *in;
	ifstream fin;
	istringstream ssin(ctfile);
	fin.open(ctfile);
	if(fin.is_open())
		in = &fin;
	else
		in = &ssin;
	
	// Start to read CT file
	j = 0;
	*in >> ct->numofbases;
	ct->allocate(ct->numofbases);
	for (ct->numofstructures = 1;!in->eof();(ct->numofstructures)++)	
	{
		ct->checknumberofstructures();
		// Read header
		strcpy (header,"");
		in->getline(header,ctheaderlength-1);
		strcpy((ct->ctlabel[ct->numofstructures]),header);
		// Read content
		for (count=1;count<=((ct->numofbases));count++)	{


			*in >> temp;//ignore base number in ctfile
			*in >> base;//read the base

			ct->nucs[count]=base[0];
			if (base[0]=='A'||base[0]=='a') ct->numseq[count]=1;
			else if (base[0]=='C'||base[0]=='c') ct->numseq[count]=2;
			else if (base[0]=='G'||base[0]=='g') ct->numseq[count]=3;
			else if (base[0]=='U'||base[0]=='u'||base[0]=='T'||base[0]=='t') ct->numseq[count]=4;
			else if (base[0]=='I') ct->numseq[count]=5;
			else ct->numseq[count]=0;
			if (ct->numseq[count]==5&&ct->numofstructures==1) {
      				ct->intermolecular = true;
       				ct->inter[j] = count;
				j++;
			}
			*in >> temp;//ignore numbering
			*in >> temp;//ignore numbering
			*in >> ct->basepr[ct->numofstructures][count];//read base pairing info
			*in >> ct->hnumber[count];//read historical numbering
			//if (ct->stacking) *in >> ct->basepr[ct->numofstructures][count+ct->numofbases];
		}
		*in >> ct->numofbases; //start on next structure and see whether the end of file is reached
	}
	(ct->numofstructures)--;
	if (fin.is_open()) fin.close();
	return;
}


//takes a nucleotide input and stores a numeric representation
void tonum(char *base,structure *ct,int count)	{
if (!strcmp(base,"A")) (ct->numseq[count] = 1);
else if(!strcmp(base,"B")) {
	(ct->numseq[count] = 1);

}
else if(!strcmp(base,"a")) {
	ct->numseq[count]=1;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"C")) (ct->numseq[count] = 2);
else if(!strcmp(base,"Z")) {
	(ct->numseq[count] = 2);

}
else if(!strcmp(base,"c")) {
	ct->numseq[count] = 2;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}
else if(!strcmp(base,"G")) (ct->numseq[count] = 3);
else if(!strcmp(base,"H")) {
	(ct->numseq[count] = 3);

}
else if(!strcmp(base,"g")) {
	ct->numseq[count] = 3;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"U")||!strcmp(base,"T")) (ct->numseq[count] = 4);
else if(!strcmp(base,"V")||!strcmp(base,"W")) {
	(ct->numseq[count] = 4);

}
else if(!strcmp(base,"u")||!strcmp(base,"t")) {
	ct->numseq[count] = 4;
	ct->nnopair++;
	ct->nopair[ct->nnopair] = count;
}

else if(!strcmp(base,"I")) {
	ct->numseq[count] = 5;
	ct->intermolecular= true;
}

else (ct->numseq[count]=0);  //this is for others, like X
return;
}

short int tonumi(char *base)	{
short int	a;
if (!strcmp(base,"A")||!strcmp(base,"B")) (a = 1);
else if(!strcmp(base,"C")||!strcmp(base,"Z")) (a = 2);
else if(!strcmp(base,"G")||!strcmp(base,"H")) (a = 3);
else if(!strcmp(base,"U")||!strcmp(base,"V")) (a = 4);
else if(!strcmp(base,"T")||!strcmp(base,"W")) (a = 4);
else (a=0);  //this is for others, like X
return a;
}


//Open seq was originally designed for reading .seq files.
//It has now been extended to automatically identify and read FASTA files,
//	where the identity line starts with a ">".

//returns 0 on error or 1 without error.

#define seqline 20000 //maximum line length in .seq file

int openseq (structure *ct,const char *seqfile) {
	char base[2],test[2],first;
	char temp[seqline],seq[seqline];
	int i,j,length,nucs;
	ifstream *FASTA;

	ct->nnopair = 0;

	//nucs = 0;

	FILE *se;

	//First check the file exists:
	//check that all the files exist with a C i/o function
	if ((se = fopen(seqfile, "r"))
		== NULL) {
		return 0;
	}



	//Now identify the file type.
	//A starting ; means .seq format.
	//A starting > means FASTA.

	//Get the first non-whitespace:
	first = (char) fgetc(se);
	while(((int) first)<33) {
		//check for premature end:
		if (first==EOF) {
			fclose(se);
			return 0;

		}
		else {
			//read the next character
			first = (char) fgetc(se);

		}

	}

	fclose(se);

	//Now we have a non whitespace character:
	if (first==';') {
		//This is a .seq file

		//read the sequence file to get the number of nucleotides
		se=fopen(seqfile,"r");

		do {
			fgets(temp,seqline,se);
			strncpy(test,temp,1);
			strcpy(test+1,"\0");
		} while (!strcmp(test,";")||((int) test[0]<32));

		//fgets(ct->ctlabel[1],seqline,se);
		strcpy(ct->ctlabel[1],temp);



		nucs = 1;
		while (1) {
			fgets(seq,seqline,se);
			length = (int) strlen (seq);
			for (j=0;j<length;j++) {//up to length-1 bcs of /n character
				strncpy (base,seq+j,1);
				strcpy (base+1,"\0");
				if (!strcmp(base,"1")) break;
				if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      				||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
					 ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
					 ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
					 !strcmp(base,"\n")||!strcmp(base,"N"))) {

						 fclose(se);
						 return 0;

				}

				if (strcmp(base," ")&&strcmp(base,"\n")) nucs++;
			}
			if (!strcmp(base,"1")) break;
		}
		ct->numofbases = nucs - 1;
		nucs--;
		fclose (se);

		if (nucs==0) return 0;

		ct->allocate(nucs);


		//now read the file
		se=fopen(seqfile,"r");

		do {
			fgets(temp,seqline,se);
			strncpy(test,temp,1);
			strcpy(test+1,"\0");
		} while (!strcmp(test,";")||((int) test[0]<32));

		//fgets(ct->ctlabel[1],seqline,se);
		strcpy(ct->ctlabel[1],temp);

		i = 1;
		while (1) {
			fgets(seq,seqline,se);
			length = (int) strlen (seq);
			for (j=0;j<length;j++) {//up to length-1 bcs of /n character
				strncpy (base,seq+j,1);
				strcpy (base+1,"\0");
				if (!strcmp(base,"1")) break;
			  if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      			||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
				 ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
				 ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
				 !strcmp(base,"\n")||!strcmp(base,"N"))) {
					 fclose(se);
					 return 0;
				 }

			  if ((!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      			||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
				 ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
				 ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base,"N"))) {


					tonum(base,ct,(i));
					ct->nucs[i]=base[0];
					ct->hnumber[i] = i;
					i++;
			  }
			}
			if (!strcmp(base,"1")) break;
		}
		//ct->numofbases = i - 1;

		fclose (se);
	}//end, this is .seq
	else if (first=='>') {
		//This is FASTA

		//read the sequence file to get the number of nucleotides
		FASTA = new ifstream(seqfile);

		//Find the opening >
		first = (char) FASTA->get();
		while(((int) first)<33) {
			//check for premature end:
			if (FASTA->eof()) {
				FASTA->close();
				return 0;

			}
			else {
				//read the next character
				first = (char) FASTA->get();

			}

		}

		//Now read the title:
		FASTA->getline(temp,seqline);//fgets(temp,seqline,se);
		strcpy(ct->ctlabel[1],temp);

		//Now remove any other lines that start with >
		FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
		while (seq[0]=='>') {
			FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
		}




		nucs = 1;
		//count the nucs
		while (1) {

			length = (int) strlen (seq);
			for (j=0;j<length;j++) {
				strncpy (base,seq+j,1);
				strcpy (base+1,"\0");

				if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      				||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
					 ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
					 ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
					 !strcmp(base,"\n")||!strcmp(base,"N"))) {

						 FASTA->close();
						 delete FASTA;
						 return 0;
				}

				if (strcmp(base," ")&&strcmp(base,"\n")) {
					nucs++;
				}
			}
			FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
			if (FASTA->eof()) {
				break;
			}

		}
		ct->numofbases = nucs - 1;
		nucs--;
		FASTA->close();//fclose (se);
		delete FASTA;

		if (nucs==0) return 0;

		ct->allocate(nucs);


		FASTA = new ifstream(seqfile);//se=fopen(seqfile,"r");

		//Find the opening >
		first = (char) FASTA->get();//(char) fgetc(se);
		while(((int) first)<33) {
			//check for premature end:
			if (FASTA->eof()) {
				fclose(se);
				return 0;

			}
			else {
				//read the next character
				first = (char) (char) FASTA->get();//fgetc(se);

			}

		}

		//Now read the title:
		FASTA->getline(temp,seqline);////fgets(temp,seqline,se);
		strcpy(ct->ctlabel[1],temp);


		FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
		while (seq[0]=='>') {
			FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
		}


		//strcpy(ct->ctlabel[0],seq);

		i = 1;

		while (i<=nucs) {

			length = (int) strlen (seq);
			for (j=0;j<length;j++) {
				strncpy (base,seq+j,1);
				strcpy (base+1,"\0");
				//if (FASTA->eof()) break;
				if (!(!strcmp(base,"A")||!strcmp(base,"a")||!strcmp(base,"C")||!strcmp(base,"c")
      				||!strcmp(base,"G")||!strcmp(base,"g")||!strcmp(base,"T")
					 ||!strcmp(base,"t")||!strcmp(base,"U")||!strcmp(base,"u")
					 ||!strcmp(base,"X")||!strcmp(base,"x")||!strcmp(base," ")||
					 !strcmp(base,"\n")||!strcmp(base,"N"))) {
						 FASTA->close();
						 delete FASTA;
						 return 0;
				}

				if (strcmp(base," ")&&strcmp(base,"\n")) {
					tonum(base,ct,(i));
					ct->nucs[i]=base[0];
					ct->hnumber[i] = i;
					i++;

				}


			}
			FASTA->getline(seq,seqline);//fgets(seq,seqline,se);
			//if (FASTA->eof()) break;

		}

		FASTA->close();//fclose (se);
		delete FASTA;



	}
	else {
		//This is neither .seq nor FASTA
		return 0;


	}


	return 1;
}

//outputs a ct file (connection table)
//Provide a pointer to structure with the pairing info, a pointer to cstring with the filename
//	if append is true, the ct table is appended an existing file, otherwise a new file is created
//	append is false by default
//By default, the columns for indicies are only 5 characters wide, so this is a problem for sequences > 9,999 nucs.
//Now, when sequences are >9,999 nucs, columns are 6 characters wide.  (The code is written yet for sequences > 99,999 nucs.)
void ctout (structure *ct,const char *ctoutfile, bool append) {
	int count,i;//length
	char line[2*ctheaderlength],number[2*numlen];//base[2]

	FILE *ctfile;
	if (append) ctfile=fopen(ctoutfile,"a");
	else ctfile=fopen(ctoutfile,"w");

	for (count=1;count<=(ct->numofstructures);count++) {

		strcpy(line,"");
		if (ct->numofbases>9999) sprintf(line,"%6i",ct->numofbases);
		else sprintf(line,"%5i",ct->numofbases);


		if (ct->energy[count]!=0) {
   			strcat(line,"  ENERGY = ");

			//gcvt((float (ct->energy[count]))/conversionfactor,6,number);
			if (conversionfactor==10)
				sprintf(number,"%.1f",(float (ct->energy[count]))/conversionfactor);
			else if (conversionfactor==100)
				sprintf(number,"%.2f",(float (ct->energy[count]))/conversionfactor);
			else sprintf(number,"%f",(float (ct->energy[count]))/conversionfactor);

   			strcat(line,number);
   			strcat(line,"  ");
		}
		else strcat(line,"  ");
		strcat(line,ct->ctlabel[count]);

		//make sure that line ends in a newline, if not, add a newline!
		if(line[strlen(line)-1]!='\n') strcat(line,"\n");

		fputs (line,ctfile);
		for (i=1;i<ct->numofbases;i++) {
			if (ct->stacking) {
				if (ct->numofbases>9999) sprintf(line,"%6i%2c%8i%6i%6i%6i%6i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
				else sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
			}
			else {
				if (ct->numofbases>9999)sprintf(line,"%6i%2c%8i%6i%6i%6i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
				else sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
					i,ct->nucs[i],(i-1),(i+1),ct->basepr[count][i],ct->hnumber[i]);
			}
			fputs(line,ctfile);
		}


		//last nucleotide not connected--
		i = ct->numofbases;
		if (ct->stacking) {
			sprintf(line,"%5i%2c%8i%5i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i],ct->basepr[count][i+ct->numofbases]);
		}
		else {
			if (ct->numofbases>9999) sprintf(line,"%6i%2c%8i%6i%6i%6i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i]); 
			
			else sprintf(line,"%5i%2c%8i%5i%5i%5i\n",
				i,ct->nucs[i],(i-1),0,ct->basepr[count][i],ct->hnumber[i]);
		}
		fputs(line,ctfile);

	}

	fclose (ctfile);
	return;
}


char *tobase (int i)

{  //function is a look up table to convert the base
	// 	integer represention to a familiar character

	if (i==1) return "A";
	 else if (i==2) return "C";
	 else if (i==3) return "G";
	 else if (i==4) return "U";
	 else if (i==0) return "X";
    else if (i==5) return "I";
	 else return "?";   //to track down mistakes

}

void sortstructures (structure *ct) {//this function reorders the structures by
													//the efn energies
register int c;
int cur,i;
char ctheader[ctheaderlength];

for (c = 2; c<=(ct->numofstructures);c++){

	cur = c;

	while (cur>1) {
		if ((ct->energy[cur])<(ct->energy[cur-1])) {
      	swap(&ct->energy[cur],&ct->energy[cur-1]);
         //also swap the ct labels:
         strcpy(ctheader, ct->ctlabel[cur]);
         strcpy(ct->ctlabel[cur],ct->ctlabel[cur-1]);
         strcpy(ct->ctlabel[cur-1],ctheader);
         for (i=1;i<=(ct->numofbases);i++) {
         	swap(&ct->basepr[cur][i],&ct->basepr[cur-1][i]);
         }
         cur--;
   	}
		else {
		break;
		}
	}
}

}


//Swap is an overloaded function that put the value of b in a and a in b.
void swap(int *a,int *b) {
int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

void swap (short int *a,short int *b) {
short int temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

void swap(float *a,float *b) {
float temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}

void swap(double *a,double *b) {
double temp;

	temp = *a;
   *a = *b;
   *b = temp;

return;

}



//#endif
