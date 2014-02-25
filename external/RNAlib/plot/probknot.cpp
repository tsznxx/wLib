
#include "probknot.h"


//Assemble the ProbKnot structure from base pair probabilities.
int ProbKnotAssemble(pfunctionclass *v, PFPRECISION *w5, structure *ct, pfdatatable *data, bool *lfce, bool *mod, PFPRECISION scaling, forceclass *fce, int iterations, int MinHelixLength) {

	PFPRECISION **probs,*rowprob;
	int i,j,iter;

	//Indicate that there will be one structure:
	ct->numofstructures = 1;

	//Build a 2-d array for storing pair probabilities, probs, note that the higher index is addressed first...
	probs = new PFPRECISION *[ct->numofbases+1];

	//also allocate space for rowprob[i], the highest probability for pairing of nucleotide i
	rowprob = new PFPRECISION [ct->numofbases+1];

	for (i=1;i<=ct->numofbases;i++) {
		probs[i] = new PFPRECISION [i+1];
		
		//Initialize rowprob to zero
		rowprob[i] = 0.0;

		//set all pairs to zero in ct
		ct->basepr[1][i] = 0;

	}

	

	
	


	
	

	//First determine pair probabilities:
	for (i=1;i<ct->numofbases;i++) {
		for (j=i+minloop+1;j<=ct->numofbases;j++) {
		
			probs[j][i] = calculateprobability(i,j,v,w5,ct,data,lfce,mod,scaling,fce);
			
			//also accumulate the best probs for each nucleotide:
			if (probs[j][i]>rowprob[i]) rowprob[i] = probs[j][i];
			if (probs[j][i]>rowprob[j]) rowprob[j] = probs[j][i];

		}
	}

	//now assemble the structure:
	for (i=1;i<ct->numofbases;i++) {
		for (j=i+minloop+1;j<=ct->numofbases;j++) {
			
			//check all possible pairs
			//take a pair if it has the highest prob for any pair involving i or j
			if (rowprob[i]==probs[j][i]&&rowprob[j]==probs[j][i]&&probs[j][i]>0.0) {

				ct->basepr[1][i]=j;
				ct->basepr[1][j]=i;
		
			}


		}
	}

	//If multiple iterations were requested, go on to do those:
	for (iter=2;iter<=iterations;iter++) {

		//starting a new iteration, reacccumulate rowprob
		for (i=1;i<=ct->numofbases;i++) rowprob[i]=0.0;
		for (i=1;i<ct->numofbases;i++) {
			for (j=i+minloop+1;j<=ct->numofbases;j++) {
			
				if (ct->basepr[1][i]==0&&ct->basepr[1][j]==0) {
				
					//accumulate the best probs for each nucleotide not already paired:
					if (probs[j][i]>rowprob[i]) rowprob[i] = probs[j][i];
					if (probs[j][i]>rowprob[j]) rowprob[j] = probs[j][i];

				}

			}
		}

		//now add to the structure:
		for (i=1;i<ct->numofbases;i++) {
			for (j=i+minloop+1;j<=ct->numofbases;j++) {
				
				if (ct->basepr[1][i]==0&&ct->basepr[1][j]==0) {
					//check all possible pairs for nucs not already paired
					//take a pair if it has the highest prob for any pair involving i or j
					if (rowprob[i]==probs[j][i]&&rowprob[j]==probs[j][i]&&probs[j][i]>0.0) {

						ct->basepr[1][i]=j;
						ct->basepr[1][j]=i;
				
					}
				}


			}
		}




	}

	//Finally, post-process the structures to remove short helices, if specified:
	if (MinHelixLength>1) {

		
		RemoveShortHelices(ct, MinHelixLength, 1);


	
	}



	//cleanup memory use:
	for (i=1;i<=ct->numofbases;i++) delete[] probs[i];
	delete[] probs;

	delete[] rowprob;

	return 0;



}



//Remove short helices, allowed stacks across single bulges
//Implemented by Stanislav Bellaousov.
void RemoveShortHelices(structure *ct, int MinHelixLength, int StructureNumber) {
	int pairs,i,j;	
	
	//Checking for helixes smaller then argv[5]
	for (i=1;i<=ct->numofbases;i++) {
	  //	  if(ct->basepr[StructureNumber][i]!=0){
	  if (ct->basepr[StructureNumber][i]>i){
	    j=ct->basepr[StructureNumber][i];
	    pairs=1;
	    while (ct->basepr[StructureNumber][i+1]==j-1||ct->basepr[StructureNumber][i+2]==j-1||ct->basepr[StructureNumber][i+1]==j-2){
	      if (ct->basepr[StructureNumber][i+1]==j-1){
			i++;
			j--;
			pairs++;
	      }
	      else if (ct->basepr[StructureNumber][i+2]==j-1){
		if (ct->basepr[StructureNumber][i+1]!=0){
		    ct->basepr[StructureNumber][ct->basepr[StructureNumber][i+1]]=0;
		    ct->basepr[StructureNumber][i+1]=0;
		  }
			i=i+2;
			j--;
			pairs++;
	      }
	      else {
			i++;
			j=j-2;
			pairs++;
	      }
	    }


	   
	    //Deleting helixes smaller then MinHelixLength

	    if (pairs<MinHelixLength){
			ct->basepr[StructureNumber][ct->basepr[StructureNumber][i]]=0;
			ct->basepr[StructureNumber][i]=0;
			
			if(i>=3){
			  while (ct->basepr[StructureNumber][i-1]==j+1||ct->basepr[StructureNumber][i-2]==j+1||ct->basepr[StructureNumber][i-1]==j+2){
			    if (ct->basepr[StructureNumber][i-1]==j+1){
			      ct->basepr[StructureNumber][ct->basepr[StructureNumber][i-1]]=0;
			      ct->basepr[StructureNumber][i-1]=0;
			      // cerr << rna->GetErrorMessage(error);
			      i--;
			      j++;
			    }
			    else if (ct->basepr[StructureNumber][i-2]==j+1){
			      ct->basepr[StructureNumber][ct->basepr[StructureNumber][i-2]]=0;
			      ct->basepr[StructureNumber][i-2]=0;
			      // cerr << rna->GetErrorMessage(error);
			      i=i-2;
			      j++;
			    }
			    else {
			      ct->basepr[StructureNumber][ct->basepr[StructureNumber][i-1]]=0;
			      ct->basepr[StructureNumber][i-1]=0;
			      // cerr << rna->GetErrorMessage(error);
			      i--;
			      j=j+2;
			      //pairs++;
			    }
			    
			  }
			}
			else if(i==2){
			  while (ct->basepr[StructureNumber][i-1]==j+1||ct->basepr[StructureNumber][i-1]==j+2){
			    if (ct->basepr[StructureNumber][i-1]==j+1){
			      ct->basepr[StructureNumber][ct->basepr[StructureNumber][i-1]]=0;
			      ct->basepr[StructureNumber][i-1]=0;
			      // cerr << rna->GetErrorMessage(error);
			      i--;
			      j++;
			    }
			    else {
			      ct->basepr[StructureNumber][ct->basepr[StructureNumber][i-1]]=0;
			      ct->basepr[StructureNumber][i-1]=0;
			      // cerr << rna->GetErrorMessage(error);
			      i--;
			      j=j+2;
			      //pairs++;
			    }
			    
			  }
			}
	    }
	  }
	}
}
