#include <cstdlib>
#include <cstring>
#include "thermodynamics.h"

//Constructor:
Thermodynamics::Thermodynamics(const bool ISRNA) {
	
	//store the backbone type
	isrna = ISRNA;
	
	//keep track of whether the parameter files have been read
	energyread=false;

	//set the default folding temperature
	temp = 310.15;

	//set the enthalpy parameters to an unread status
	enthalpy = NULL;


}

//Set the folding temperature:
//Return an error code pertaining to reading the thermodynamic parameters:
int Thermodynamics::SetTemperature(double temperature) {
	//int error;

	//set the folding temperature
	temp = temperature;

	//If the thermodynamic parameter files were read at some point, delete them now:
	if (energyread) delete data;

	//Setting energyread to false will ensure that the parameters will be re-read from disk
		//and set for the correct temperature at ay point they are needed.
	energyread = false;

	return 0;

	//now read the thermodynamic parameters:
	//error = ReadThermodynamic();
	//if (error!=0) return error;
	//else {
	//	return 0;
	//}



}

//return the current folding temperature
double Thermodynamics::GetTemperature() {
	return temp;
}



//Destructor:
Thermodynamics::~Thermodynamics() {

	//If the thermodynamic parameter files were read at some point, delete them now:
	if (energyread) delete data;

	//If the enthalpy parameters were read from disk, they must be deleted now:
	if (enthalpy!=NULL) delete enthalpy;

}

/*	Function GetDat

	Function gets the names of data files to open

*/

void Thermodynamics::GetDat(char *loop, char *stackf, char *tstackh, char *tstacki,
	char *tloop, char *miscloop, char *danglef, char *int22,
	char *int21,char *coax, char *tstackcoax,
    char *coaxstack, char *tstack, char *tstackm, char *triloop,
    char *int11, char *hexaloop, char *tstacki23, char *tstacki1n,
	char *datapath, bool isRNA, bool isEnthalpy)

{

 
  //copy the path to each name
  strcpy (loop,datapath);
  strcpy (stackf,datapath);
  strcpy (tstackh,datapath);
  strcpy (tstacki,datapath);
  strcpy (tloop,datapath);
  strcpy (miscloop,datapath);
  strcpy (danglef,datapath);
  strcpy (int22,datapath);
  strcpy (int21,datapath);
  strcpy (triloop,datapath);
  strcpy (coax,datapath);
  strcpy (tstackcoax,datapath);
  strcpy (coaxstack,datapath);
  strcpy (tstack,datapath);
  strcpy (tstackm,datapath);
  strcpy (int11,datapath);
  strcpy (hexaloop,datapath);
  strcpy (tstacki23,datapath);
  strcpy (tstacki1n,datapath);


  if( !isRNA) {
	  //these are dna parameters and so they need to start with "dna"
	strcat (loop,"dna");
	strcat (stackf,"dna");
	  strcat (tstackh,"dna");
	  strcat (tstacki,"dna");
	  strcat (tloop,"dna");
	  strcat (miscloop,"dna");
	  strcat (danglef,"dna");
	  strcat (int22,"dna");
	  strcat (int21,"dna");
	  strcat (triloop,"dna");
	  strcat (coax,"dna");
	  strcat (tstackcoax,"dna");
	  strcat (coaxstack,"dna");
	  strcat (tstack,"dna");
	  strcat (tstackm,"dna");
	  strcat (int11,"dna");
	  strcat (hexaloop,"dna");
	  strcat (tstacki23,"dna");
	  strcat (tstacki1n,"dna");  
	 
  }
	 
	//append the actual file name
  strcat (loop,"loop.");
  strcat (stackf,"stack.");
  strcat (tstackh,"tstackh.");
  strcat (tstacki,"tstacki.");
  strcat (tloop,"tloop.");
  strcat (miscloop,"miscloop.");
  strcat (danglef,"dangle.");
  strcat (int22,"int22.");
  strcat (int21,"int21.");
  strcat (triloop,"triloop.");
  strcat (coax,"coaxial.");
  strcat (tstackcoax,"tstackcoax.");
  strcat (coaxstack,"coaxstack.");
  strcat (tstack,"tstack.");
  strcat (tstackm,"tstackm.");
  strcat (int11,"int11.");
  strcat (hexaloop,"hexaloop.");
  strcat (tstacki23,"tstacki23.");
  strcat (tstacki1n,"tstacki1n.");
  
  if (isEnthalpy) {
	  //thse are enthalpy partameters so they need to end in .dh
	strcat (loop,"dh");
	strcat (stackf,"dh");
	  strcat (tstackh,"dh");
	  strcat (tstacki,"dh");
	  strcat (tloop,"dh");
	  strcat (miscloop,"dh");
	  strcat (danglef,"dh");
	  strcat (int22,"dh");
	  strcat (int21,"dh");
	  strcat (triloop,"dh");
	  strcat (coax,"dh");
	  strcat (tstackcoax,"dh");
	  strcat (coaxstack,"dh");
	  strcat (tstack,"dh");
	  strcat (tstackm,"dh");
	  strcat (int11,"dh");
	  strcat (hexaloop,"dh");
	  strcat (tstacki23,"dh");
	  strcat (tstacki1n,"dh");  
  }
  else {
	  //these are free energy parameters and the files end in .dat
	strcat (loop,"dat");
	strcat (stackf,"dat");
	  strcat (tstackh,"dat");
	  strcat (tstacki,"dat");
	  strcat (tloop,"dat");
	  strcat (miscloop,"dat");
	  strcat (danglef,"dat");
	  strcat (int22,"dat");
	  strcat (int21,"dat");
	  strcat (triloop,"dat");
	  strcat (coax,"dat");
	  strcat (tstackcoax,"dat");
	  strcat (coaxstack,"dat");
	  strcat (tstack,"dat");
	  strcat (tstackm,"dat");
	  strcat (int11,"dat");
	  strcat (hexaloop,"dat");
	  strcat (tstacki23,"dat");
	  strcat (tstacki1n,"dat");  


  }


}

datatable *Thermodynamics::GetDatatable() {

	return data;

}



//read the thermodynamic parameters from disk at location $DATAPATH or pwd
int Thermodynamics::ReadThermodynamic(const char *pathname) {
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;

	//only allocate the datatable if energyread is false, meaning that no parameters are loaded
	//	This is important because the user might alter the temperature with SetTemperature(), triggering a re-read of the parameters.
	if (!energyread) data = new datatable();
	

	//Set the path to the thermodynamic parameters:
	if (pathname!=NULL) {
		//The user is specifying a path to the thermodynamic parameters
		strcpy(datapath,pathname);
		strcat(datapath,"/");

	}
	else {
		//Get the path to thermodynamic parameters from $DATAPATH, if available
		pointer = getenv("DATAPATH");
		if (pointer!=NULL) {
			strcpy(datapath,pointer);
			strcat(datapath,"/");
		}
		else strcpy(datapath,"");

	}
	
	//open the data files -- must reside in pwd or $DATAPATH.
	//open the thermodynamic data tables
	GetDat (loop, stackf, tstackh, tstacki,tloop, miscloop, danglef, int22,
          int21,coax, tstackcoax,coaxstack, tstack, tstackm, triloop,
          int11, hexaloop, tstacki23, tstacki1n, datapath, isrna);//the true indicates RNA parameters
	if (opendat (loop,stackf,tstackh,tstacki,tloop,miscloop,danglef,int22,int21,
   		coax,tstackcoax,coaxstack,tstack,tstackm,triloop,int11,hexaloop,tstacki23, tstacki1n,data)==0) {
      	

		delete data;
		energyread = false;//energy files have not been correctly read
		return 5;//an error code
		
	}
	else {
		
		//now check to see if the temperature is other than 310.15 K:
		if (temp>(310.15+TOLERANCE)||temp<(310.15-TOLERANCE)) {

			//temperature is altered, so the enthalpy tables need to be read:
			datatable *localenthalpy;

			//get the names of the enthalpy files
			GetDat(loop, stackf, tstackh, tstacki,
				tloop, miscloop, danglef, int22,
				int21,coax, tstackcoax,
				coaxstack, tstack, tstackm, triloop,
				int11, hexaloop, tstacki23, tstacki1n, datapath, isrna,true);//rtue means this is enthalpy files

			//allocate a table to storte enthalpy parameters
			localenthalpy = new datatable();

			//open the enthlpy parameters and check for errors
			if (opendat(loop, stackf, tstackh, tstacki,
				tloop, miscloop, danglef, int22,
				int21,coax, tstackcoax,
				coaxstack, tstack, tstackm, triloop,
				int11,hexaloop,tstacki23, tstacki1n, localenthalpy)==0) {
	
				delete localenthalpy;
				return 5;//an error has occured

			}

			//using the enthalpy parameters and the folding free energy changes at 37 degrees C,
			//	calculate folding free energy changes at temp, the current temperature, and overwrite the 
			//	folding free energy changes at 37 degrees C
			dG_T((float)temp,*data,*localenthalpy,*data);

			//remove the enthalpy parameters, there is no need to keep them
			delete localenthalpy;

		}


		energyread=true;//energy files have been correctly read
		return 0;
	}
}


// This function is used to provide an enthalpy table.
datatable *Thermodynamics::GetEnthalpyTable() {
	char loop[maxfil],stackf[maxfil],tstackh[maxfil],tstacki[maxfil],
		tloop[maxfil],miscloop[maxfil],danglef[maxfil],int22[maxfil],
		int21[maxfil],coax[maxfil],tstackcoax[maxfil],
		coaxstack[maxfil],tstack[maxfil],tstackm[maxfil],triloop[maxfil],int11[maxfil],hexaloop[maxfil],
		tstacki23[maxfil], tstacki1n[maxfil],datapath[maxfil],*pointer;

	//start by determining if the parameters need to be read or whether they have already been read:
	if (enthalpy==NULL) {
		//The parameters have not been read, so read them now:

		//Get the information from $DATAPATH, if available
		pointer = getenv("DATAPATH");
		if (pointer!=NULL) {
			strcpy(datapath,pointer);
			strcat(datapath,"/");
		}
		else strcpy(datapath,"");

		//get the names of the enthalpy files
		GetDat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11, hexaloop, tstacki23, tstacki1n, datapath, isrna,true);//true means this is enthalpy files

		//allocate a table to storte enthalpy parameters
		enthalpy = new datatable();

		//open the enthlpy parameters and check for errors
		if (opendat(loop, stackf, tstackh, tstacki,
			tloop, miscloop, danglef, int22,
			int21,coax, tstackcoax,
			coaxstack, tstack, tstackm, triloop,
			int11,hexaloop,tstacki23, tstacki1n, enthalpy)==0) {

			delete enthalpy;
			return NULL;//an error has occured

		}
		else return enthalpy;

		

	}
	else {
		//The parameters have been read, so return the pointer to enthalpy
		return enthalpy;

	}


}

//Copy thermodynamic parameters from an instance of Thermodynamics class.

void Thermodynamics::CopyThermodynamic(Thermodynamics *thermo) {

	if (thermo->GetEnergyRead()) {
		SetTemperature(thermo->GetTemperature());
		data = new datatable();
		energyread = true;
		

		*data = *thermo->GetDatatable(); 
	}

	return;

}


//Return whether this instance of Thermodynamics has the paremters populated (either from disk or from another Thermodynamics class).
		
		
bool Thermodynamics::GetEnergyRead() {

	return energyread;

}

