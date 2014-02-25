
//Class RNA -- Wrapper for RNAstructure for use in object-oriented applications

#include "RNA.h"

#include "structure.h"
#include "arrayclass.h"
#include "dotarray.h"
#include "forceclass.h"
#include "rna_library.h"
#include "stackclass.h"
#include "stackstruct.h"



#include "algorithm.h"
#include "outputconstraints.h"
#include "pfunction.h"
#include "boltzmann.h"
#include "alltrace.h"
#include "stochastic.h"
#include "MaxExpect.h"
#include "probknot.h"

#include <iostream>

const float epsilon = 1e-6; // a small number for a tolerance in comparing floats


//constructor where user provides a string with the sequence
RNA::RNA(const char sequence[], const bool IsRNA):Thermodynamics(IsRNA) {
	int i,j;



	//allocate ct
	ct = new structure();

	//Specify the sequence length based on the string length
	ct->numofbases = (short) strlen(sequence);
	ct->allocate(ct->numofbases);//allocate the space required in the arrays for the sequence


	//Now store the sequence information
	for (i=1;i<=ct->numofbases;i++) {
		if (sequence[i-1]=='A'||sequence[i-1]=='a') ct->numseq[i]=1;
		else if (sequence[i-1]=='C'||sequence[i-1]=='c') ct->numseq[i]=2;
		else if (sequence[i-1]=='G'||sequence[i-1]=='g') ct->numseq[i]=3;
		else if (sequence[i-1]=='U'||sequence[i-1]=='u'||sequence[i-1]=='T'||sequence[i-1]=='t') ct->numseq[i]=4;
		else ct->numseq[i]=0;

		ct->nucs[i] = sequence[i-1];
		ct->hnumber[i]=i;
	}

	//Also make it safe for the user to not specify all the structural information,e.g. set all the pairs to zero:
	for (i=1;i<ct->allocatedstructures;i++) {
		for (j=1;j<=ct->numofbases;j++) {

			ct->basepr[i][j]=0;
		}

	}

	//no structures have been specified
	ct->numofstructures=0;

	//no lables have been specified
	strcpy(ct->ctlabel[1],"\n");


	//These should not be located here: (DHM 5/3/2011)
	// Experimental bonuses array
   	//ct->EX = new double *[ct->numofbases+1];
   	//for(i=0;i<ct->numofbases+1;i++) {
	//  ct->EX[i] = new double[ct->numofbases+1];
   	//}
	//for (i=0;i<ct->numofbases+1;i++){
	//  for (j=0;j<ct->numofbases+1;j++){
	 //   ct->EX[i][j] = 0.0;
	 // }
     //   }


	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//set error status to zero
	ErrorCode=0;

	//Do not report progress by default:
	progress=NULL;



}

//constructor where user provides a string with a filename
//	The flag type indicates the file type: type = 1 => ct file, type = 2 => .seq file, type = 3 => .pfs file.
//	The fconstructor saves an error code in ErrorCode:
//	0 = no error, 1 = file not found
//	2 = error opening file.
//  GetErrorCode provides public access to this errorcode.
RNA::RNA(const char filename[], const int type, const bool IsRNA ):Thermodynamics(IsRNA) {


	//allocate ct
	ct = new structure();

	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not (yet) read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//Do not report progress by default:
	progress=NULL;

	ErrorCode = FileReader(filename, type);
	return;
}

//Default constructor.
RNA::RNA(const bool IsRNA):Thermodynamics(IsRNA) {





	//allocate the underlying structure class and nothing more.
	//User is then required to propogate the sequence and structural information manually.
	ct = new structure();

	//set the number of nucleotides to zero because no sequence has been read
	ct->numofbases = 0;

	//Indicate that the partition function calculation has not been performed.
	partitionfunctionallocated = false;

	//Indicate that the energy data is not (yet) read.
	energyallocated = false;

	//Drawing coordinates have not been determined.
	drawallocated = false;

	//set error status to zero
	ErrorCode=0;

	//Do not report progress by default:
	progress=NULL;



}


//Return the value of ErrorCode
int RNA::GetErrorCode() {
	return ErrorCode;
}

//Return a c string that describes errors from GetErrorCode and other errors.
char* RNA::GetErrorMessage(const int error) {

	if (error==0) return "No Error.\n";
	else if (error==1) return "Input file not found.\n";
	else if (error==2) return "Error opening file.\n";
	else if (error==3) return "Structure number out of range.\n";
	else if (error==4) return "Nucleotide number out of range.\n";
	else if (error==5) return "Error reading thermodynamic parameters.\nPlease set environment variable $DATAPATH to the location of the thermodynamic parameters.\n";
	else if (error==6) return "This would form a pseudoknot and is not allowed.\n";
	else if (error==7) return "This pair is non-canonical and is therefore not allowed.\n";
	else if (error==8) return "Too many restraints specified.\n";
	else if (error==9) return "This nucleotide already under a conflicting constraint.\n";
	else if (error==10) return "There are no structures to write to file.\n";
	else if (error==11) return "Nucleotide is not a U.\n";
	else if (error==12) return "Maximum pairing distance is too short.\n";
	else if (error==13) return "Error reading constraint file.\n";
	else if (error==14) return "A traceback error occurred.\n";
	else if (error==15) return "No partition function data is available.\n";
	else if (error==16) return "Wrong save file version used or file format not recognized.\n";
	else if (error==17) return "This function cannot be performed unless a save file (.sav) was correctly loaded by the RNA constructor.\n";
	else if (error==18) return "This threshold is too low to generate valide secondary structures.\n";
	else if (error==19) return "The structure coordinates have not been determined, use DetermineDrawingCoordinates() to calculate the coordinates.\n";
	else if (error==20) return "No sequence has been read.\n";
	else if (error==21) return "Probabilities summed to greater than 1 in stochastic traceback.\n";
	else if (error==22) return "Programming error.  Incorrect file type passed to constructor.\n";
	else if (error==23) return "There are no structures present.\n";
	else if (error==24) return "Too few iterations.  There must be at least one iteration.\n";
	else if (error==25) return "Index is not a multiple of 10.\n";
	else if (error==26) return "k, the equilibrium constant, needs to be greater than or equal to 0.\n";
	else return "Unknown Error\n";


}

//Return a string that describes errors from GetErrorCode and other errors.
//This uses GetErrorMessage to acrually get the errors.
std::string RNA::GetErrorMessageString(const int error) {
	std::string temp;

	temp = GetErrorMessage(error);

	return temp;


}


//User specifies a base pair between i and j in structure # structurenumber.
int RNA::SpecifyPair(const int i, const int j, const int structurenumber) {
	int a,b;


	//start with error checking:
	if (i<0||i>ct->numofbases||j<0||j>ct->numofbases) return 4;
	else if (structurenumber<1) return 3;

	//also keep track of the maximum number of structures for which pairs have been specified
	if (structurenumber>ct->numofstructures) {
		if (structurenumber>=ct->allocatedstructures) {
			//Expand the space allocated for structures
			int **tempbasepr;
			char **tempctlabel;
			int *tempenergy;
			int tempallocatedstructures;

			//allocate pointers for tracking currently used memory
			tempbasepr = new int *[ct->allocatedstructures];
			tempctlabel = new char *[ct->allocatedstructures];

			//store pointers to currently used memory
			tempenergy = ct->energy;
			for (a=0;a<ct->allocatedstructures;a++) {
				tempbasepr[a] = ct->basepr[a];
				tempctlabel[a] = ct->ctlabel[a];
			}

			tempallocatedstructures = ct->allocatedstructures;

			//set the number of structures to 1000 greater than the current number
			ct->allocatedstructures = structurenumber+1000;

			//have new memory alllocated
			ct->allocatestructure();

			//now copy over all the information from the old memory
			for (a=1;a<=ct->numofstructures;a++) {
				ct->energy[a] = tempenergy[a];
				strcpy(ct->ctlabel[a],tempctlabel[a]);
				for (b=1;b<=ct->numofstructures;b++) {
					ct->basepr[a][b]=tempbasepr[a][b];


				}
				//delete old memory
				delete[] tempctlabel[a];
				delete[] tempbasepr[a];

			}
			delete[] tempenergy;
			delete[] tempctlabel;
			delete[] tempbasepr;

		}
		for (a=ct->numofstructures+1;a<=structurenumber;a++) {

			//also set the comment to be empty
			strcpy(ct->ctlabel[a],"\n");


			for (b=1;b<=ct->numofbases;b++) {
				ct->basepr[a][b]=0;
			}

		}


		ct->numofstructures = structurenumber;

	}

	//now register the pair:
	ct->basepr[structurenumber][i]=j;
	ct->basepr[structurenumber][j]=i;


	return 0;

}
//Break a pair that i is involved in
int RNA::RemoveBasePair(const int i, const int structurenumber) {

	//start with error checking:
	if (i<0||i>ct->numofbases) return 4;
	else if (structurenumber<1||structurenumber>=ct->allocatedstructures) return 3;


	//break the pair for i and i's pairing partner
	if (ct->basepr[structurenumber][i]!=0) {


		ct->basepr[structurenumber][ct->basepr[structurenumber][i]]=0;
		ct->basepr[structurenumber][i]=0;
	}

	//return that there was no error
	return 0;

}

//remove all pairs in structure # structurenumber.
//Also, roll back the number of specified structures if this is the last specified structure.
int RNA::RemovePairs(const int structurenumber) {
	int i;

	//do some error checking
	if (structurenumber>ct->numofstructures||structurenumber<1) return 5;

	//reset pairs
	for (i=1;i<=ct->numofbases;i++) ct->basepr[structurenumber][i]=0;

	//decrement the number of structures, if appropriate
	if (structurenumber==ct->numofstructures)ct-> numofstructures--;

	return 0;

}


//Calculate and return the folding free energy change for structure number structurenumber.
double RNA::CalculateFreeEnergy(const int structurenumber, const bool UseSimpleMBLoopRules) {
	//Do some simple error checking
	if (structurenumber<1||structurenumber>ct->numofstructures) return 0.0;

	if (!energyread) {
		//The thermodynamic data tables have not yet been read
		if (ReadThermodynamic()!=0) {
			ErrorCode = 5;//record an error
			return 0.0;//return 0.0 if a problem occurs

		}
		else ErrorCode=0;//record that there was no error.
	}
	else ErrorCode=0;//Set the error code to zero because no errors were encountered.


	efn2(data,ct,structurenumber,UseSimpleMBLoopRules);

	//conversion factor is set in defines.h.  Free energies are multiplied by this factor internally so that integer math can be used.
	return (((double) ct->energy[structurenumber])/conversionfactor);


}

//Write the details on the energy caclulation for all structures.
int RNA::WriteThermodynamicDetails(const char filename[], const bool UseSimpleMBLoopRules) {

	if (!energyread) {
		//The thermodynamic data tables have not yet been read
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs
	}
	efn2(data,ct,0,UseSimpleMBLoopRules,filename);
	return 0;


}


//Predict the secondary structure by free energy minimization.
//Also generate subooptimal solutions using a heuristic.
int RNA::FoldSingleStrand(const float percent, const int maximumstructures, const int window, const char savefile[], const int maxinternalloopsize, bool mfeonly) {
	char *savefilename;
	int percenti;
	int tracebackstatus;

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;


	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}

	//savefile will be passed to the function dynamic for structure prediction.
	//Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
	if (!strcmp(savefile,"")) savefilename = NULL;
	else {
		savefilename=new char[((int) strlen(savefile))+1];
		strcpy(savefilename,savefile);
	}

	//right now, dynamic requires an integer specification of percent energy change.
	//FoldSingleStrand takes this is a float to provide the opportunity to reform this in the future.
	//For now, cast percent as an integer to pass to dynamic.
	percenti = (int) percent;

	//Predict the secondary structures.
	tracebackstatus = dynamic(ct, data, maximumstructures, percenti, window, progress, false, savefilename, maxinternalloopsize,mfeonly);

	//Clean up the memory use.
	delete[] savefilename;

	if(tracebackstatus!=0) return 14;//This indicates a traceback error.
	else return 0;



}

// Predict the lowest free energy secondary structure and generate all suboptimal structures.
int RNA::GenerateAllSuboptimalStructures(const float percent, const double deltaG) {

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;


	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}


	//Call the alltrace function to do the work:
	alltrace(ct,data, ((short) percent), ((short) (deltaG*conversionfactor)),progress,NULL);

	return 0;


}

// Predict the structure with maximum expected accuracy and suboptimal structures.
int RNA::MaximizeExpectedAccuracy(const double maxPercent, const int maxStructures, const int window, const double gamma) {

	//first trap some possible errors
	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping
	MaxExpectFill(ct, v, w5, pfdata, lfce, mod, fce, maxPercent, maxStructures, window, gamma, progress);

	return 0;//no error return functionality right now

}

// This function predicts structures composed of probable base pairs.
int RNA::PredictProbablePairs(const float probability) {
	int i,j,count;
	char thresh[8];

	//first trap some possible errors
	if (probability > epsilon && probability < 0.500-epsilon) {
		//The threshold is too low to be valie and not low enough that it will be considered zero, a the default
		return 18;
	}

	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping

	//First reset the pairs to zero:
	for (i=1;i<=ct->numofbases;i++) {

			ct->basepr[1][i] = 0;

	}
	if (probability>epsilon) {
		//The user specified a threshold, so use that and generate one structure
		ct->numofstructures = 1;



		for (i=1;i<ct->numofbases;i++) {
			for (j=i+1;j<=ct->numofbases;j++) {

				if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce) > probability) {
					//This pair exceeded the threshold, so add it to the list
					ct->basepr[1][i] = j;
					ct->basepr[1][j] = i;

				}

			}
		}

		//put the threshold in the ctlabel, which will appear in the header of a ct file
		sprintf(thresh,"%1.5f",probability);
		strcpy(ct->ctlabel[2],ct->ctlabel[1]);
		strcpy(ct->ctlabel[1]," >");
		strcat(ct->ctlabel[1],thresh);
		strcat(ct->ctlabel[1]," pairing probability; ");
		strcat(ct->ctlabel[1],ct->ctlabel[2]);



	}
	else {
		//The default threshold was specified, so create 8 structures, with thresholds of >=0.99, >=0.97, >=0.95, >=0.90, >=0.80, >=0.70, >=0.60, >0.50.

		ct->numofstructures = 8;

		//zero the basepair arrays for structures 2 through 8
		for (i=1;i<=ct->numofbases;i++) {

				ct->basepr[2][i] = 0;
				ct->basepr[3][i] = 0;
				ct->basepr[4][i] = 0;
				ct->basepr[5][i] = 0;
				ct->basepr[6][i] = 0;
				ct->basepr[7][i] = 0;
				ct->basepr[8][i] = 0;
		}

		//loop over the structures and thresholds
		for (count=1;count<=8;count++) {

			for (i=1;i<ct->numofbases;i++) {
				for (j=i+1;j<=ct->numofbases;j++) {

					if (count==1) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.99) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==2) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.97) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==3) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.95) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==4) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.90) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==5) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.80) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==6) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.70) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==7) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>=.60) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
					else if (count==8) {
						if (calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce)>.50) {
							ct->basepr[count][i]=j;
							ct->basepr[count][j]=i;
						}
					}
				}
			}

		}

		//add lables that would appear in a ct file header
		strcpy(ct->ctlabel[9],ct->ctlabel[1]);
		strcpy(ct->ctlabel[2]," >=97% probable pairs ");
		strcpy(ct->ctlabel[3]," >=95% probable pairs ");
		strcpy(ct->ctlabel[4]," >=90% probable pairs ");
		strcpy(ct->ctlabel[5]," >=80% probable pairs ");
		strcpy(ct->ctlabel[6]," >=70% probable pairs ");
		strcpy(ct->ctlabel[7]," >=60% probable pairs ");
		strcpy(ct->ctlabel[8]," >50% probable pairs ");
		strcat(ct->ctlabel[2],ct->ctlabel[1]);
		strcat(ct->ctlabel[3],ct->ctlabel[1]);
		strcat(ct->ctlabel[4],ct->ctlabel[1]);
		strcat(ct->ctlabel[5],ct->ctlabel[1]);
		strcat(ct->ctlabel[6],ct->ctlabel[1]);
		strcat(ct->ctlabel[7],ct->ctlabel[1]);
		strcat(ct->ctlabel[8],ct->ctlabel[1]);
		strcpy(ct->ctlabel[1]," >=99% probable pairs ");
		strcat(ct->ctlabel[1],ct->ctlabel[9]);

	}

	return 0;



}

//Calculate the partition function for the current sequence.
int RNA::PartitionFunction(const char savefile[], double temperature) {
	int i,j;
	char *savefilename;

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	if (!energyread) {
		//The thermodynamic data tables have not been read and need to be read now.
		if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

	}

	//savefile will be passed to the function dynamic for structure prediction.
	//Dynamic expects a null pointer if no file is to be created, so set savefilename to null if savefile is an empty string.
	if (!strcmp(savefile,"")) savefilename = NULL;
	else {
		savefilename=new char[((int) strlen(savefile))+1];
		strcpy(savefilename,savefile);
	}

	if (partitionfunctionallocated) {
		delete v;
		delete w;
		delete wmb;
		delete wl;
		delete wmbl;
		delete wcoax;
		delete fce;
		delete lfce;
		delete mod;
		delete w3;
		delete w5;
		delete pfdata;



	}

	//Allocate the memory needed (only if this is the first call to pfunction):

	//indicate that the memory has been allocated so that the destructor will delete it.
	partitionfunctionallocated = true;

	//allocate space for the v and w arrays:
	w = new pfunctionclass(ct->numofbases);
	v = new pfunctionclass(ct->numofbases);
	wmb = new pfunctionclass(ct->numofbases);
	wl = new pfunctionclass(ct->numofbases);
	wmbl = new pfunctionclass(ct->numofbases);
	wcoax = new pfunctionclass(ct->numofbases);
	fce = new forceclass(ct->numofbases);

	lfce = new bool [2*ct->numofbases+1];
	mod = new bool [2*ct->numofbases+1];



	for (i=0;i<=2*ct->numofbases;i++) {
		lfce[i] = false;
		mod[i] = false;
	}


	for (i=1;i<=ct->nmod;i++) {

		if (ct->mod[i]!=1&&ct->mod[i]!=ct->numofbases) {
			mod[ct->mod[i]]=true;
			mod[ct->mod[i]+ct->numofbases]=true;
		}
	}




	w5 = new PFPRECISION [ct->numofbases+1];
	w3 = new PFPRECISION [ct->numofbases+2];


	if (ct->intermolecular) {
		//take advantage of templating to prevent intramolecular base pairs

		ct->allocatetem();//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
		for (i=1;i<ct->inter[0];i++) {
			for (j=i+1;j<=ct->inter[2];j++) {

				//Set intermolecular pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
				ct->tem[j][i]=false;

			}
		}
		for (i=ct->inter[2]+1;i<ct->numofbases;i++) {
			for (j=i+1;j<=ct->numofbases;j++) {

				//Another set of intermolecular pairs to forbid.
				ct->tem[j][i]=false;

			}
		}



	}

	//Initialize the partition function datatable:
		//Ignore the setting of parameter temperature if it is less than zero.
		//Generally, this parameter should be left at the default.
	if (temperature<0) pfdata=new pfdatatable(data,scalingdefinition,temp);
	else pfdata=new pfdatatable(data,scalingdefinition,temperature);

	//This code converts the SHAPE array of data to equilibrium constants.  This is
	//needed for the partition function.  NOTE, however, there is no going back so a structure
	//that has been used for partition functions cannot then be used to predict a structure
	//by free energy minimization. This is a compromise for efficiency, but clearly something
	//that could cause a problem for programmers.
	if (ct->shaped) {
	  for (i=1;i<=2*ct->numofbases;i++) ct->SHAPE[i]=boltzman( ct->SHAPE[i], pfdata->temp);
	}

	if (ct->experimentalPairBonusExists){
	  //symmetrize as well...
	  for (i = 1; i <= 2*ct->numofbases; i++){
	    for(j = i; j <= 2*ct->numofbases; j++){
	      double avg = (double)( 0.5*(ct->EX[i][j] + ct->EX[j][i]) );
	      ct->EX[i][j] = boltzman( avg, pfdata->temp);
	      ct->EX[j][i] = ct->EX[i][j];
	    }
	  }
	}


	//The next section handles the case where base pairs are not
			//not allowed to form between nucs more distant
			//than ct->maxdistance
	if (ct->limitdistance) {

		//This allocates a bool array that indicates what pairs are allowed.  It is initilaized to true, i.e. all nucs can pair with each other.
		if (!ct->templated) ct->allocatetem();

		for (j=minloop+2;j<=ct->numofbases;j++) {
			for (i=1;i<j;i++) {
				//Set distant pairs to false, i.e. not allowed.  Note the indexing with the high index and then the low index number.
				if (j-i>=ct->maxdistance) ct->tem[j][i]=false;
			}
		}



	}


	//This is the workhorse function:
	calculatepfunction(ct,pfdata,progress,savefilename,false,&Q,w,v,wmb,wl,wmbl,wcoax,fce,w5,w3,mod,lfce);


	if (savefilename!=NULL) {
		writepfsave(savefilename,ct,w5,w3,v,w,wmb,wl,wmbl,wcoax,fce,mod,lfce,pfdata);

		//clean up some memory use:
		delete[] savefilename;
	}



	return 0;

}

//Predict maximum expected accuracy structures that contain pseudoknots.


int RNA::ProbKnot(int iterations, int MinHelixLength) {

	//first trap some possible errors
	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}

	if (iterations < 1) {
		//there can't be fewer than one iteration
		return 24;

	}



	//Past error trapping
	//Call the ProbKnot Program:
	return ProbKnotAssemble(v, w5, ct, pfdata, lfce, mod, pfdata->scaling, fce, iterations, MinHelixLength );


}


//Refold a sequence using data from a save file.
int RNA::ReFoldSingleStrand(const float percent, const int maximumstructures, const int window) {

	if (!energyallocated) {
		//A .sav file was not read by the constructor.  Therefore, this function cannot be performed.
		return 17;

	}

	//Now do the refolding.
	return traceback(ct, data, ev, ew, ewmb, ew2, ewmb2, ew3, ew5, fce, lfce, vmin, maximumstructures, (int) percent, window,mod);

}


//Sample structures from the Boltzman ensemble.
int RNA::Stochastic(const int structures, const int seed) {

	if (!partitionfunctionallocated) {
		//There is no partition function data available.
		return 15;
	}


	//Past error trapping, call the stochastic traceback function
	return stochastictraceback(w,wmb,wmbl,wcoax,wl,v,
		fce, w3,w5,pfdata->scaling, lfce, mod, pfdata, structures,
		ct, seed, progress);


}


//Force a nucleotide to be double stranded (base paired).
//Return an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceDoubleStranded(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->numofbases) return 4;//i is out of range

	//Check for conflicting restraints (anything specifying i to be unpaired):
	for (index=1;index<=ct->nnopair;index++) {

		if(ct->nopair[index]==i)	return 9;
	}

	if (ct->ndbl>=maxforce) {
		return 8;//Too many restraints specified.
	}
	else {
		(ct->ndbl)++;
		ct->dbl[ct->ndbl] = i;
		return 0;

	}


}

//Function to specify a nucleotide, i, that is a U in a GU pair
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint, 11 = nucleotide not U).
int RNA::ForceFMNCleavage(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->numofbases) return 4;//i is out of range

	//Check to make sure the nucleotide is U.
	if (ct->numseq[i]!=4) return 11;

	//Check for conflicting restraints (anything specifying i to be unpaired):
	for (index=1;index<=ct->nnopair;index++) {

		if(ct->nopair[index]==i) return 9;
	}

	//Check for a nucleotide already forced to be in a pair that is not a GU pair.
	for (index=1;index<=ct->npair;index++) {

		if (i==ct->pair[index][0]&&ct->numseq[ct->pair[index][1]]!=3) return 9;
		else if (i==ct->pair[index][1]&&ct->numseq[ct->pair[index][0]]!=3) return 9;

	}


	if (ct->ngu>=maxforce) {
		return 8;//Too many restraints specified.
	}
	else {
		(ct->ngu)++;
		ct->gu[ct->ngu] = i;
		return 0;

	}


}

//Specify the maximum distance allowed between paired nucleotides in subsequent structure prediction.
//return An integer that indicates an error code (0 = no error, 12 = too long or too short).
int RNA::ForceMaximumPairingDistance(const int distance) {

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	if (distance < minloop+1) return 12;
	else {
		ct->limitdistance = true;
		ct->maxdistance=distance;
		return 0;
	}


}


//Indicate a nucleotide that is accessible to chemical modification.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified).
int RNA::ForceModification(const int i) {

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->numofbases) return 4;//i is out of range


	//check if too many modifications have been specified.
	if (ct->nmod>=maxforce) return 8;



	else {
		//Go ahead and record the constraint.
		ct->nmod++;
		ct->mod[ct->nmod]=i;
		return 0;
	}


}

//Force a base pair between nucleotides i and j.
//Returns an error code: (0 = no error, 4 = nucleotide out of range, 6 = pseudoknot formation, 7 = non-canonical pair, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForcePair(const int i, const int j) {
	bool allowedpairs[6][6]={{false,false,false,false,false,false},{false,false,false,false,true,false},{false,false,false,true,false,false},{false,false,true,false,true,false},
	{false,true,false,true,false,false},{false,false,false,false,false,false}};
	int index;
	int locali,localj;

	//First perform the error checking:

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//Note: In structure, forced pairs run between index of 1 and a maximum of maxforce-1.

	if (ct->npair==(maxforce-1)) return 8;//This means there are too many pair constraints.

	if (i<1||i>ct->numofbases) return 4;//i is out of range
	if (j<1||j>ct->numofbases) return 4;//j is out of range

	if (!allowedpairs[ct->numseq[i]][ct->numseq[j]]) return 7;//non-canonical pair

	//sort indexes from 5' to 3':
	if (i>j) {
		locali=j;
		localj=i;
	}
	else {
		locali=i;
		localj=j;
	}

	//check for pseudoknots with any other forced pair or the same nucleotide forced into two pairs:
	for (index=1;index<=ct->npair;index++) {
		if (locali<ct->pair[index][0]&&ct->pair[index][0]<localj&&localj<ct->pair[index][1]) return 6;//a pseudoknot

		if (locali==ct->pair[index][0]||locali==ct->pair[index][1]||localj==ct->pair[index][0]||localj==ct->pair[index][1]) return 9;//i or j is in another forced pair

	}

	//now check for other conflicting restraints:
	for (index=0;index<ct->nforbid;index++) {
		if(ct->forbid[index][0]==locali && ct->forbid[index][1]==localj ) return 9;//The pair was forbidden.

	}
	for (index=1;index<=ct->nnopair;index++) {

		if(ct->nopair[index]==locali||ct->nopair[index]==localj)	return 9;//i or j was previously forced single-stranded.
	}

	//Now register the restraint because the error checking was clear or errors.
	ct->npair++;
	ct->pair[ct->npair][0]=locali;
	ct->pair[ct->npair][1]=localj;

	return 0;

}

//Prohibit a pair between two nucleotides in subsequent structure prediction.
//Returns an integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = nucleotide in conflicting restraint).
int RNA::ForceProhibitPair(const int i, const int j) {
	int index,locali,localj;

	//First perform the error checking:

	//Note: In structure, forced pairs run between index of 0 and a maximum of maxforce-1.

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	if (ct->nforbid==(maxforce)) return 8;//This means there are too many pair constraints.

	if (i<1||i>ct->numofbases) return 4;//i is out of range
	if (j<1||j>ct->numofbases) return 4;//j is out of range


	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}

	//check to make sure this pair hasn't been forced:
	for (index=1;index<=ct->npair;index++) {

		if (locali==ct->pair[index][0]&&localj==ct->pair[index][1]) return 9;//i or j is in a forced pair

	}


	//Now register the restraint because the error checking was clear or errors.

	ct->forbid[ct->nforbid][0]=locali;
	ct->forbid[ct->nforbid][1]=localj;
	ct->nforbid++;

	return 0;

}

//Force a nucleotide to be single stranded in subsequent structure prediction.
//An integer that indicates an error code (0 = no error, 4 = nucleotide out of range, 8 = too many restraints specified, 9 = same nucleotide in conflicting restraint).
int RNA::ForceSingleStranded(const int i) {
	int index;

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//Check that nucleotide is valid
	if (i<1||i>ct->numofbases) return 4;//i is out of range

	//Check for conflicting constraints; anything forcing a nucleotide to be paired.

	for (index=1;index<=ct->npair;index++) {//check all the forced pairs
		if (i==ct->pair[index][0]||i==ct->pair[index][1]) return 9;//i is in a forced pair
	}
	for (index=1;index<=ct->ndbl;index++) {//check all the force doubles
		if(ct->dbl[index]==i)	return 9;
	}
	for (index=1;index<=ct->ngu;index++) {//check all the force FMN
		if(ct->gu[index]==i)	return 9;
	}


	//Check if too many restraints have been specified
	if (ct->nnopair>=maxforce) return 8;
	else {

		ct->nnopair++;
		ct->nopair[ct->nnopair]=i;

	}

	return 0;

}

//Return a nucleotide that is forced double stranded.
int RNA::GetForcedDoubleStranded(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->ndbl) return 0;

	//Now return the constraint.
	return ct->dbl[constraintnumber+1];//note that the underlying ct indexes from 1 to ndbl.

}

//Return a nucleotide that is accessible to FMN cleavage.
int RNA::GetForcedFMNCleavage(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->ngu) return 0;

	//Now return the constraint.
	return ct->gu[constraintnumber+1];//note that the underlying ct indexes from 1 to ndbl.

}

//Return a nucleotide that is accessible to modification.
int RNA::GetForcedModification(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->nmod) return 0;

	//Now return the constraint.
	return ct->mod[constraintnumber+1];//note that the underlying ct indexes from 1 to ndbl.

}

//Return a nucleotide in a forced pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedPair(const int constraintnumber, const bool fiveprime) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->npair) return 0;

	//Now return the constraint.
	if (fiveprime) return ct->pair[constraintnumber+1][0];//note that the underlying ct indexes from 1 to ndbl.
	else return ct->pair[constraintnumber+1][1];

}

//Return a nucleotide in a prohibited pair.
//fiveprime determines if the nucleotide is the five prime or the three prime nucleotide in the constraint.  true = five prime nucleotide.
int RNA::GetForcedProhibitedPair(const int constraintnumber, const bool fiveprime) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->nforbid) return 0;

	//Now return the constraint.
	if (fiveprime) return ct->forbid[constraintnumber][0];//note that the underlying ct indexes from 0 to ndbl-1.
	else return ct->forbid[constraintnumber][1];

}

//Return a nucleotide that is forced single stranded.
int RNA::GetForcedSingleStranded(const int constraintnumber) {

	//First make sure the constraintnumber is valid.
	if (constraintnumber<0||constraintnumber>ct->nnopair) return 0;

	//Now return the constraint.
	return ct->nopair[constraintnumber+1];//note that the underlying ct indexes from 1 to ndbl.


}


//Return the maximum pairing distance.
//Return an integer that indicates the maximum distance allowed between paired nucleotides, where -1 indicates that the maximum distance is not set.
int RNA::GetMaximumPairingDistance() {

	if (ct->limitdistance) return ct->maxdistance;
	else return -1;

}


//Return the number of nucletides forced to be paired.
int RNA::GetNumberOfForcedDoubleStranded() {

	return ct->ndbl;

}

// Add an experimental bonus to a pair of nucleotides
//void RNA::SetExperimentalBonus(const int i, const int j, const double bonus){

//	ct->EX[i][j] = bonus;

//}


//!Return the number of nucleotides accessible to FMN cleavage.
int RNA::GetNumberOfForcedFMNCleavages() {

	return ct->ngu;

}

//!Return the number of nucleotides accessible to chemical modification.
int RNA::GetNumberOfForcedModifications() {

	return ct->nmod;

}

//!Return the number of forced base pairs.
int RNA::GetNumberOfForcedPairs() {

	return ct->npair;

}

//!Return the number of prohibited base pairs.
int RNA::GetNumberOfForcedProhibitedPairs() {

	return ct->nforbid;

}

//!Return the number of nucleotides that are not allowed to pair.
int RNA::GetNumberOfForcedSingleStranded() {

	return ct->nnopair;

}

//Read a set of folding constraints to disk in a plain text file.
//filename is a c string that is the file name to be read.
//Returns an integer that indicates an error code (0 = no error, 1 = file not found, 13 = error reading constraint file).
int RNA::ReadConstraints(const char filename[]) {
	FILE *check;

	//check that the file exists.
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}
	fclose(check);

	//Now read the constraints
	if (readconstraints(filename, ct)) return 0;
	else return 13;



}

//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//IsPseudoEnergy indicates whether this is the pseudo folding free energy constraint (the preferred method).  This defaults to true.
//parameter1 is the slope when IsPseudoEnergy=true and is a threshold above which nucleotides are forced single stranded otherwise.
//parameter2 is the intercept when IsPseudoEnergy=true and is a threshold above which a nucleotide is considered chemically modified otherwise.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, diffSHAPE, DMS, and CMCT). Defaults to SHAPE.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, std::string modifier, const bool IsPseudoEnergy) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	fclose(check);

	if (IsPseudoEnergy) {
		//This is the pseudo energy version

		ct->SHAPEslope=parameter1*conversionfactor;//register the slope in tenths of kcal/mol
		ct->SHAPEintercept=parameter2*conversionfactor;//register the intercept in tenths of a kcal/mol
		ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies

	}
	else {
		ct->ReadSHAPE(filename,(float) parameter1,(float) parameter2);//call ReadSHAPE() with parameters to parse thresholds
	}
	return 0;


}


//Read SHAPE data to constrain structure prediction on subsequent structure predictions.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the slope.
//parameter2 is the intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadExperimentalPairBonus(const char filename[], double const experimentalOffset, double const experimentalScaling ) {
	FILE *check;

	//check that the SHAPE input file exists
	if ( strlen( filename ) > 0  ) {
	  if ( (check = fopen(filename, "r"))== NULL) {
	    //the file is not found
	    fclose(check);
	    return 1;
	  }

	  fclose(check);
	}

	ct->ReadExperimentalPairBonus(filename, experimentalOffset, experimentalScaling );

	return 0;

}


//Read SHAPE data to constrain structure prediction on subsequent structure predictions - overloaded version for including single-stranded SHAPE.
//filename is a c string that indicates a file that contains SHAPE data.
//parameter1 is the double-stranded slope.
//parameter2 is the double-stranded intercept.
//modifier is the type of chemical modification probe that was used (currently accepted values are SHAPE, DMS, and CMCT). Defaults to SHAPE.
//ssm is the single-stranded slope.
//ssb in the single-stranded intercept.
//Returns an integer that indicates an error code (0 = no error, 1 = input file not found).
int RNA::ReadSHAPE(const char filename[], const double parameter1, const double parameter2, const double ssm, const double ssb, std::string modifier) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	fclose(check);



	ct->SHAPEslope=parameter1*conversionfactor;//register the slope in tenths of kcal/mol
	ct->SHAPEintercept=parameter2*conversionfactor;//register the intercept in tenths of a kcal/mol
	ct->SHAPEslope_ss=ssm*conversionfactor;//register the slope in tenths of kcal/mol
	ct->SHAPEintercept_ss=ssb*conversionfactor;//register the intercept in tenths of a kcal/mol
	ct->ReadSHAPE(filename, modifier);//call ReadSHAPE() to read the file and determine pseudo energies


	return 0;


}

//Read Double Strand Offset
int RNA::ReadDSO(const char filename[]) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	ct->ReadOffset(NULL,filename);

	return 0;

}

//Read Single Strand Offset
int RNA::ReadSSO(const char filename[]) {
	FILE *check;

	//check that the SHAPE input file exists
	if ((check = fopen(filename, "r"))== NULL) {
		//the file is not found
		fclose(check);
		return 1;
	}

	ct->ReadOffset(filename,NULL);

	return 0;

}


//Remove all previously defined constraints.
void RNA::RemoveConstraints() {

	ct->ndbl = 0;
	ct->npair = 0;
	ct->nnopair = 0;
	ct->nmod = 0;
	ct->ngu=0;
	ct->nforbid=0;
	ct->min_gu=0;
	ct->min_g_or_u=0;
	ct->nneighbors=0;
	ct->nregion=0;
	ct->nmicroarray=0;

}

//Add extrinsic restraints for partition function calculations.
int RNA::SetExtrinsic(int i, int j, double k) {
	int locali, localj;


	//First do the error checking:
	//check the indices
	if (i<1||i>ct->numofbases||j<1||j>ct->numofbases) return 4;

	//make sure the equilibrium constant is not less than zero.
	if (k<0) return 26;

	//Now past error checking:

	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}



	if (ct->constant==NULL) {
		//allocate the space needed in structure to store the constants
		ct->allocateconstant();
	}


	ct->constant[localj][locali] = k;

	return 0;
}


//Write the set of folding constraints to disk.
int RNA::WriteConstraints(const char filename[]) {

	outputconstraints(filename,ct);
	return 0;

}



//Specify a comment for inclusion in subsequently written .ct files.
int RNA::AddComment(const char comment[], const int structurenumber) {

	//start with error checking:
	if (structurenumber<1||structurenumber>=ct->allocatedstructures) return 3;

	//now register the comment (at the end of existing string, but before the newline there by default):
	strcpy(ct->ctlabel[structurenumber]+(strlen(ct->ctlabel[structurenumber])-1),comment);

	//add a newline at the end of the comment -- required when ct is written
	strcat(	ct->ctlabel[structurenumber],"\n");

	//also keep track of the maximum number of structures for which pairs have been specified
	if (structurenumber>ct->numofstructures) ct->numofstructures = structurenumber;
	return 0;

}

//Write a ct file of the structures
int RNA::WriteCt(const char filename[], bool append) {
	if (ct->numofstructures>0) {
		ctout(ct,filename,append);
		return 0;
	}
	else return 10; //an error code


}

//Write a ct file of the structures
int RNA::WriteDotBracket(const char filename[]) {
	if (ct->numofstructures>0) {
		ct->writedotbracket(filename);
		return 0;
	}
	else return 10; //an error code


}


//Break any pseudoknots that might be in a structure.
int RNA::BreakPseudoknot(const bool minimum_energy, const int structurenumber) {

	int i,j,structures;
	structure *tempct;


	double **bpProbArray; //contains the raw bp probabilities for each bp
	double *bpSSProbArray; //contains the raw single strand probability for a base
	double **vwArray;  //contains v and w recursion values
	double **vwPArray; //the v' and w' recursion values

	double *w3Array=0;//w3Array[i] is the maximum score from nucletides i to ct->numofbases
	double *w5Array=0;//w5Array[i] is the maximum score from nucleotides 1 to i
							

	int Length;
	int start,stop;


	double sumPij; //holds the sum of probabilities for base pairs based on a specific i over js


	//Make sure there are structures...
	if (ct->numofstructures<=0) {
		return 23;
	}

	//If a specific structure is requested, i.e. structurenumber!=0, make sure the structure is valid.
	if (structurenumber!=0) {
		if (structurenumber<0||structurenumber>ct->numofstructures) return 3;

	}

	if (minimum_energy) {
		//mnimum_energy==true, so minimize free energy in the pseudoknot free structure


		if (!energyread) {
			//The thermodynamic data tables have not been read and need to be read now.
			if (ReadThermodynamic()!=0) return 5;//return non-zero if a problem occurs

		}

		//allocate tempct:
		tempct = new structure(2);
		tempct->allocate(ct->numofbases);
		strcpy(tempct->ctlabel[1],"temp\n");
		tempct->numofbases = ct->numofbases;

		for (i=1;i<=ct->numofbases;i++) {
			tempct->numseq[i]=ct->numseq[i];

		}
		//allocate a template of allowed pairs
		tempct->allocatetem();

		//Break pairs for each structure (first check if only one structure should be treated)

		if (structurenumber!=0) {
			start = structurenumber;
			stop = structurenumber;

		}
		else {
			start = 1;
			stop = ct->numofstructures;
		}
		for (structures=start;structures<=stop;structures++) {

			//initialize all base pairs as unallowed:
			for (i=0;i<=ct->numofbases;i++) {
				for (j=i+1;j<=ct->numofbases;j++) {
    				tempct->tem[j][i] = false;
				}
		   }

			//now allow the pairs that are in the loaded ct:
			for (i=1;i<=ct->numofbases;i++) {
				if (ct->basepr[structures][i]>i) {
					tempct->tem[ct->basepr[structures][i]][i] = true;
				}
			}

			tempct->numofstructures = 0; //strip the structure from the ct at this point

			//Predict the secondary structures.
			alltrace(tempct,data,0,0,progress,NULL,false);//
			//dynamic(tempct, data, 1, 0, 0, progress, false, NULL, 30);



			//copy the pairs back to ct
			for (i=1;i<=ct->numofbases;i++) {
				ct->basepr[structures][i]=tempct->basepr[1][i];

			}

			//Also copy the energy back to ct
			ct->energy[structures]=tempct->energy[1];



		}
		delete tempct;

	}//end of minimum_energy == true
	else {
		//This is minimum_energy == false, so maximize pairs

		//allocate tempct:
		tempct = new structure(2);
		tempct->allocate(ct->numofbases);
		strcpy(tempct->ctlabel[1],"temp\n");
		tempct->numofbases = ct->numofbases;

		for (i=1;i<=ct->numofbases;i++) {
			tempct->numseq[i]=ct->numseq[i];

		}


		//loop over all structures unless the programmer requested a specific structure
		if (structurenumber!=0) {
			start = structurenumber;
			stop = structurenumber;

		}
		else {
			start = 1;
			stop = ct->numofstructures;
		}
		for (structures=start;structures<=stop;structures++) {



			tempct->numofstructures = 0; //strip the structure from the ct at this point




			//allocate main arrays and initialize the Arrays to defaults
			bpProbArray = new double *[ct->numofbases+1];
			bpSSProbArray = new double [ct->numofbases+1];
			vwArray = new double *[ct->numofbases+1];
			vwPArray = new double *[ct->numofbases+1];


			sumPij = 0;

			for (i=0;i<=ct->numofbases;i++) {
				bpProbArray[i] = new double [ct->numofbases+1];
				vwArray[i] = new double [ct->numofbases+1];
				vwPArray[i] = new double [ct->numofbases+1];

				tempct->basepr[1][i] = 0;
				bpSSProbArray[i] = 0;

				for (j=0;j<=ct->numofbases;j++) {
					bpProbArray[i][j]=0;
					vwArray[i][j]=-0;
					vwPArray[i][j]=-0;
				}
			}




			tempct->nucs[0] = ' ';

			// Recursion rules investigate for vwArray
			//    1)  if the base pair (BP) probability value is 0, skip that pair
			//    2)  hairpin turns at 5 BP
			//    3)  stack/internal/bulge pairing at 7 BPs
			//    4)  multibranching at 12 BPs (2 hairpins and a stack)
			// Because of the rules for hairpin, the probabilities
			//    can be taken from the bpProbArray directly until j-i > 5

			// Calculate the single stranded probabilities for each base
			// Pi = 1 - (for all j, sum(Pij)
			// fill in w for the diagonal for the Pi,i
			for (i=1; i<=ct->numofbases; i++)
			{


				if (ct->basepr[structures][i]!=0) bpSSProbArray[i] = 0.0;
				else bpSSProbArray[i] = 1.0;



				vwArray[i][i] = bpSSProbArray[i];
			} // end loop over each base pair


			//Calculate the base pair probabilities to start...
			for (Length = 2; Length <=ct->numofbases; Length++)
			{

				//begin populating v and w along the diagonal starting with the
				//   shortest hairpin length
				for (i=1, j=i+Length-1; j<=ct->numofbases; i++, j++)
				{
					if (ct->basepr[structures][i]==j) {
						bpProbArray[j][i]=1.0;
					}
					else {
						bpProbArray[j][i]=-1.0;
					}

				}
			}

			//Call the MEAFill routine.
				//Note the false at the end "allows" non-canonical pairs.  This is required so that
				//non-canonical pairs aren't spuriosly broken
			MEAFill(tempct, bpProbArray, bpSSProbArray, vwArray, vwPArray, w5Array, w3Array, 1.0, 0, progress,false);



			// start traceback
			trace(tempct, vwArray, vwPArray, bpProbArray, 1.0, 0, 1, 0);






			// Deallocate memory for the MaxExpect calculation
			//Arrays with functionality in the fill step
			for (i=0;i<=ct->numofbases;i++) {
				delete[] bpProbArray[i];
			}
			delete[] bpProbArray;

			delete[] bpSSProbArray;

			for (i=0; i<=ct->numofbases; i++) {
				delete[] vwArray[i];
				delete[] vwPArray[i];
			}
			delete[] vwArray;
			delete[] vwPArray;




			//copy the pairs back to ct
			for (i=1;i<=ct->numofbases;i++) {
				ct->basepr[structures][i]=tempct->basepr[1][i];

			}

			//Also copy the energy back to ct
			ct->energy[structures]=tempct->energy[1];



		}
		delete tempct;



	}
	return 0;





}


// Report if there are any pseudoknots in a structure.


bool RNA::ContainsPseudoknot(const int structurenumber) {
	int i,j;

	//make sure structurenumber is a valid structure
	if (structurenumber<1||structurenumber>ct->numofstructures) {
		ErrorCode=3;
		return false;
	}
	else {
		//passed error trapping:

		//check all nucs for pairs
		for (i=1;i<ct->numofbases;i++) {

			if (ct->basepr[structurenumber][i]>i) {
				//found pair, check for crossing pairs
				for (j=i+1;j<ct->basepr[structurenumber][i];j++) {

					if (ct->basepr[structurenumber][j]>j) {

						if (ct->basepr[structurenumber][j]>ct->basepr[structurenumber][i]) {

							return true;

						}

					}
				}
			}
		}

		return false;//no pseudoknot was found
	}
}


//Get the ensemble folding free energy change as determined by the partition function.
double RNA::GetEnsembleEnergy() {

	//check to see if partition function data is present.
	if (!partitionfunctionallocated) {
		ErrorCode = 15;
		return 0.0;
	}

	//past error trapping, so set no error found
	ErrorCode = 0;

	//calculate the ensemble folding free energy change, -RT ln (Q).
	//One needs to also account for the fact that a scaling is applied per nucleotide.
	return (double) ((-RKC*temp)*(log(w5[ct->numofbases])-ct->numofbases*log(pfdata->scaling)));

}

//Get the folding free energy change for a predicted structure.
double RNA::GetFreeEnergy(const int structurenumber) {

	//make sure structurenumber is a valid structure
	if (structurenumber<1||structurenumber>ct->numofstructures) {
		ErrorCode=3;
		return 0.0;
	}
	else {
		//error trapping complete
		return (((double) ct->energy[structurenumber])/((double) conversionfactor));

	}

}

// Get the nucleotide to which the specified nucleotide is paired.
int RNA::GetPair(const int i, const int structurenumber) {

	//check to make sure i is a valid nucleotide
	if (i<1||i>ct->numofbases) {
		ErrorCode=4;
		return 0;
	}
	//now make sure structurenumber is a valid structure
	else if (structurenumber<1||structurenumber>ct->numofstructures) {
		ErrorCode=3;
		return 0;
	}
	else {
		//error trapping complete
		return ct->basepr[structurenumber][i];

	}

}

//Extract the lowest free energy for a structure that contains the i-j pair using data from a save file (.sav).
double RNA::GetPairEnergy(const int i, const int j) {
	int locali, localj;

	//check whether the a save file was read (by the constructor)
	if (!energyallocated) {
		ErrorCode = 17;
		return 0.0;
	}


	if (i<1||i>ct->numofbases) {
		//i is out of range
		ErrorCode = 4;
		return 0.0;
	}
	if (j<1||j>ct->numofbases) {
		//j is out of range
		ErrorCode = 4;
		return 0.0;
	}


	//sort indexes from 5' to 3':
	if (i>j) {
		locali = j;
		localj = i;
	}
	else {
		locali=i;
		localj=j;
	}

	//No error;
	ErrorCode = 0;

	//calculate and return the energy
	return ((((double) (ev->f(locali,localj)+ev->f(localj,locali+ct->numofbases)))/conversionfactor));



}

//Get the total number of specified structures
int RNA::GetStructureNumber() {

		return ct->numofstructures;

}

//return the base pairing probability between nucleotides i and j
double RNA::GetPairProbability(const int i, const int j) {

	//check to see if partition function data is present.
	if (!partitionfunctionallocated) {
		ErrorCode = 15;
		return 0.0;
	}

	//check that the nucleotides are in the correct range
	if (i<1||j>ct->numofbases||j<0||j>ct->numofbases) {
		ErrorCode = 4;
		return 0.0;
	}

	//past error trapping, so set no error found
	ErrorCode = 0;

	//calculate the base pair probability
	return (double) calculateprobability(i,j,v,w5,ct,pfdata,lfce,mod,pfdata->scaling,fce);


}

//Determine the coordinates for drawing a secondary structure.
int RNA::DetermineDrawingCoordinates(const int height, const int width, const int structurenumber) {

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) return 20;

	//First check for errors in the specification of the structure number:
	if (structurenumber<0||structurenumber>ct->numofstructures) return 3;

	if (!drawallocated) {
		//If the memory has not been allocated for coordinates, go ahead and allocate it now
		structurecoordinates = new coordinates(ct->numofbases);
		drawallocated = true;
	}

	//now perform the calculation for coordinate determination:
	place(structurenumber, ct, structurecoordinates, height, width);

	//no errors, return 0
	return 0;


}

//Provide the comment from the ct file as a string.
std::string RNA::GetCommentString(const int structurenumber) {
	std::string temp;


	//start with error checking:
	if (structurenumber<1||structurenumber>=ct->allocatedstructures){
		//The request is for a structure that is out of range
		ErrorCode = 3;
		temp = "";

		return temp;

	}

	temp = ct->ctlabel[structurenumber];

	return temp;

}

//Get the X coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideXCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->numofbases) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->x[i];


}

//Get the Y coordinate for nucleotide i for drawing a structure.
int RNA::GetNucleotideYCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->numofbases) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->y[i];


}

// Get the X coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelXCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->numofbases) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	if (i%10!=0) {
		//The nucleotide index is not a multiple of 10, return an error
		ErrorCode = 25;
		return 0;

	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->num[i/10][0];


}

// Get the Y coordinate for placing the nucleotide index label specified by i.
int RNA::GetLabelYCoordinate(const int i) {

	if (!drawallocated) {
		//The drawing coordinates were not pre-determined by DetermineDrawingCoordinates(), so indicate an error.
		ErrorCode = 19;
		return 0;
	}

	if (i<0||i>ct->numofbases) {
		//The nucleotide is invalid, indicate an error.
		ErrorCode = 4;
		return 0;
	}

	if (i%10!=0) {
		//The nucleotide index is not a multiple of 10, return an error
		ErrorCode = 25;
		return 0;

	}

	//Fetch the coordinate from the coordinates structure, structurecoordinates.
	return structurecoordinates->num[i/10][1];

}

//Get the identity of nucleotide i.
char RNA::GetNucleotide(const int i) {

	//check to make sure that a sequence has been read
	if (ct->numofbases==0) {
		ErrorCode = 20;
		return '-';
	}

	//Check that nucleotide is valid
	if (i<1||i>ct->numofbases) {
		ErrorCode = 4;//i is out of range
		return '-';
	}

	return ct->nucs[i];

}

//Get the total length of the sequence
int RNA::GetSequenceLength() {

	return ct->numofbases;

}

//Return the type of backbone (true = RNA, false = DNA).
bool RNA::GetBackboneType() {

	return isrna;

}






//Access the underlying structure class.
structure *RNA::GetStructure() {

	return ct;

}


//Provide a TProgressDialog for following calculation progress.
//A TProgressDialog class has a public function void update(int percent) that indicates the progress of a long calculation.
void RNA::SetProgress(TProgressDialog& Progress) {


	progress = &Progress;

	return;

}


//Provide a means to stop using a TProgressDialog.
//StopProgress tells the RNA class to no longer follow progress.  This should be called if the TProgressDialog is deleted, so that this class does not make reference to it.
void RNA::StopProgress() {

	progress=NULL;
	return;

}

TProgressDialog* RNA::GetProgress() {

		return progress;

}


RNA::~RNA() {






	if (partitionfunctionallocated) {
		//The partition function calculation was performed, so the memory allocated for the partition function needs to be deleted.
		delete[] lfce;
		delete[] mod;
		delete[] w5;
		delete[] w3;
		delete v;
		delete w;
		delete wmb;
		delete wl;
		delete wmbl;
		delete wcoax;
		delete fce;
		delete pfdata;

	}

	if (energyallocated) {
		//A folding save file was opened, so clean up the memory use.

		delete[] lfce;
		delete[] mod;



		delete[] ew5;
		delete[] ew3;


		if (ct->intermolecular) {
			delete ew2;
			delete ewmb2;
		}

		delete ev;
		delete ew;
		delete ewmb;
		delete fce;

	}

	if (drawallocated) {
		//The drawing coordiantes have been determined, so the memory needs to be cleaned up.
		delete structurecoordinates;

	}

	delete ct;//delete the structure


}



//This is a protected function for handling file input.
int RNA::FileReader(const char filename[], const int type) {
	FILE *check;
	short vers;
	//double scaling;
	int i;

	//open the file based on type:
	if (type==1) {
		//type indicates a ct file
		openct(ct,filename);//TODO: add error handling for ct opening
		return 0;

	}
	else if (type==2) {
		//type indicates a .seq file

		ct->numofstructures=0;
		if (openseq(ct,filename)==1) return 0;//no error
		else return 2;//File open error



	}
	else if (type==3) {
		//type indicates a partition function save file

		//allocate the ct file by reading the save file to get the sequence length:
		std::ifstream sav(filename,std::ios::binary);

		read(&sav,&(vers));//read the version of the save file

		if (vers!=pfsaveversion) {
			//Wrong version!
			sav.close();
			return 16;
		}

		//read the length of the sequence
		read(&sav,&(ct->numofbases));
		sav.close();

		//allocate everything
		ct->allocate(ct->numofbases);

		w = new pfunctionclass(ct->numofbases);
		v = new pfunctionclass(ct->numofbases);
		wmb = new pfunctionclass(ct->numofbases);
		wmbl = new pfunctionclass(ct->numofbases);
		wcoax = new pfunctionclass(ct->numofbases);
		wl = new pfunctionclass(ct->numofbases);
		fce = new forceclass(ct->numofbases);

		w5 = new PFPRECISION [ct->numofbases+1];
		w3 = new PFPRECISION [ct->numofbases+2];

		lfce = new bool [2*ct->numofbases+1];
		mod = new bool [2*ct->numofbases+1];

		pfdata = new pfdatatable();

		//indicate that the memory has been allocated so that the destructor will delete it.
		partitionfunctionallocated = true;

		//load all the data from the pfsavefile:
		readpfsave(filename, ct, w5, w3,v, w, wmb,wl, wmbl, wcoax, fce,&pfdata->scaling,mod,lfce,pfdata);
		return 0;

	}
	else if (type==4) {
		//Type indicates restoration of a folding save file.


		//peek at the length of the sequence and whether the folding is intermolecular to allocate arrays:
		std::ifstream sav(filename,std::ios::binary);

		//Read the save file version information.
		read(&sav,&vers);

		//check the version
		if (vers!=safiversion) {
			//Wrong version!
			sav.close();
			return 16;
		}

		read(&sav,&(ct->numofbases));
		read(&sav,&(ct->intermolecular));
		sav.close();

		//indicate that everything is allocated and needs to be deleted in the destructor
		energyallocated = true;

		//allocate everything
		ct->allocate(ct->numofbases);

		ew = new arrayclass(ct->numofbases);
		ev = new arrayclass(ct->numofbases);
		ewmb = new arrayclass(ct->numofbases);
		fce = new forceclass(ct->numofbases);


		lfce = new bool [2*ct->numofbases+1];
		mod = new bool [2*ct->numofbases+1];

		ew5 = new integersize [ct->numofbases+1];
		ew3 = new integersize [ct->numofbases+2];

		if (ct->intermolecular) {
			ew2 = new arrayclass(ct->numofbases);
			ewmb2 = new arrayclass(ct->numofbases);

			for (i=0;i<3;i++) read(&sav,&(ct->inter[i]));

		}
		else {
			ew2 = NULL;
			ewmb2 = NULL;
		}

		//indicate that the thermodynamic parameters are read and available (and need to be deleted).
		energyread = true;
		data = new datatable();

		//now read the file.
		readsav(filename, ct, ew2, ewmb2, ew5, ew3, lfce, mod, data,
				 ev, ew, ewmb, fce, &vmin);

		//set error status to zero
		 return 0;


	}
	return 22;
}




//The following should not be included for compilations for Windows:
#ifndef _WINDOWS_GUI

//A global function for error reporting
void errmsg(int err,int erri) {

if (err==30) {
	std::cout << "End Reached at traceback #"<<erri<<"\n";
   return;
}
if (err==100) {
	std::cout << "error # "<<erri;
   return;
}
switch (err) {
	case 1:
   	std::cout << "Could not allocate enough memory";
      break;
   case 2:
   	std::cout << "Too many possible base pairs";
      break;
   case 3:
   	std::cout << "Too many helixes in multibranch loop";
   case 4:
   	std::cout << "Too many structures in CT file";
   default:
   	std::cout << "Unknown error";
}
//std::cin >> err;
return;

}

#endif
