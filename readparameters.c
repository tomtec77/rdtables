/**
 * @file readparameters.c
 * @brief Read input parameters from file
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"

char Path1[256];

double BaryonFrac;

#ifdef TESTMODE
double Omega0, OmegaLambda;
double UnitMass_in_g, UnitLength_in_cm, UnitVelocity_in_cm_per_s;
double UnitMass_in_Msun, UnitVelocity_in_km_per_s, UnitLength_in_kpc;
#endif

int NumOutputTimes;

double TimeOfFinalOutput;
double TimeBetOutputs;

int NumCbins;

int NumMdiscBins;

double StartingMdisc;
double WidthOfMdiscBins;

int NumLambdaBins;

double StartingLambda;
double WidthOfLambdaBins;

int NumMbulgeBins;

double StartingMbulge;
double WidthOfMbulgeBins;

/*double MassFractionForVc;*/
double OpticalRadius;

unsigned int OutputListOn;
unsigned int ConcentrationListOn;
unsigned int MbulgeListOn;
unsigned int BulgeAngularMomentumOn;
unsigned int RdTablesOn;
unsigned int VcTablesOn;

char NameOutputList[256];
char NameConcentrationList[256];
char NameMbulgeList[256];
char NameConcentrationTables[256];


void readparameterfile(char *filename)
{
  FILE *fp;
  char buffer[256], buffer2[256];
  char *error;
  int elements;
  
  fp = fopen(filename, "r");
  if (fp==NULL) {
    fprintf(stdout, "Error: cannot open file '%s' for reading: %s (%u)\n", 
	    filename, strerror(errno), errno);
    exit(EXIT_FAILURE);     
  }
  fprintf(stdout,"Reading from parameters file '%s'...\n", filename);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Path1%s", Path1);
  checkforerror(error,elements,&buffer[0]);
  printf("    Path1                       %s\n", Path1);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BaryonFrac%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  BaryonFrac = atof(buffer2);
  printf("    BaryonFrac                  %g\n", BaryonFrac);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NumOutputTimes%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  NumOutputTimes = atoi(buffer2);
  printf("    NumOutputTimes              %d\n", NumOutputTimes);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeOfFinalOutput%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  TimeOfFinalOutput = atof(buffer2);
  printf("    TimeOfFinalOutput           %g\n", TimeOfFinalOutput);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"TimeBetOutputs%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  TimeBetOutputs = atof(buffer2);
  printf("    TimeBetOutputs              %g\n", TimeBetOutputs);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NumCbins%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  NumCbins = atoi(buffer2);
  printf("    NumCbins                    %d\n", NumCbins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NumMdiscBins%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  NumMdiscBins = atoi(buffer2);
  printf("    NumMdiscBins                %d\n", NumMdiscBins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"StartingMdisc%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  StartingMdisc = atof(buffer2);
  printf("    StartingMdisc               %g\n", StartingMdisc);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"WidthOfMdiscBins%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  WidthOfMdiscBins = atof(buffer2);
  printf("    WidthOfMdiscBins            %g\n", WidthOfMdiscBins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NumLambdaBins%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  NumLambdaBins = atoi(buffer2);
  printf("    NumLambdaBins               %d\n", NumLambdaBins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"StartingLambda%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  StartingLambda = atof(buffer2);
  printf("    StartingLambda              %g\n", StartingLambda);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"WidthOfLambdaBins%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  WidthOfLambdaBins = atof(buffer2);
  printf("    WidthOfLambdaBins           %g\n", WidthOfLambdaBins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NumMbulgeBins%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  NumMbulgeBins = atoi(buffer2);
  printf("    NumMbulgeBins               %d\n", NumCbins);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"StartingMbulge%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  StartingMbulge = atof(buffer2);
  printf("    StartingMbulge              %g\n", StartingMbulge);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"WidthOfMbulgeBins%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  WidthOfMbulgeBins = atof(buffer2);
  printf("    WidthOfMbulgeBins           %g\n", WidthOfMbulgeBins);

  /*error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MassFractionForVc%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  MassFractionForVc = atof(buffer2);
  printf("    MassFractionForVc           %g\n", MassFractionForVc);*/
  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"OpticalRadius%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);   
  OpticalRadius = atof(buffer2);
  printf("    OpticalRadius               %g\n", OpticalRadius);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"OutputListOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  OutputListOn = atoi(buffer2);
  printf("    OutputListOn                %d\n", OutputListOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"ConcentrationListOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  ConcentrationListOn = atoi(buffer2);
  printf("    ConcentrationListOn         %d\n", ConcentrationListOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"MbulgeListOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  MbulgeListOn = atoi(buffer2);
  printf("    MbulgeListOn                %d\n", MbulgeListOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"BulgeAngularMomentumOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  BulgeAngularMomentumOn = atoi(buffer2);
  printf("    BulgeAngularMomentumOn      %d\n", BulgeAngularMomentumOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"RdTablesOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  RdTablesOn = atoi(buffer2);
  printf("    RdTablesOn                  %d\n", RdTablesOn);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"VcTablesOn%s",buffer2);  
  checkforerror(error,elements,&buffer[0]);
  VcTablesOn = atoi(buffer2);
  printf("    VcTablesOn                  %d\n", VcTablesOn);

  printf("\n");

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameOutputList%s",NameOutputList);
  checkforerror(error,elements,&buffer[0]);
  printf("    NameOutputList              %s\n", NameOutputList);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameMbulgeList%s",NameMbulgeList);
  checkforerror(error,elements,&buffer[0]);
  printf("    NameMbulgeList              %s\n", NameMbulgeList);

  /* Tables are cosmology-independent, but to produce values in physical
     units when running test mode a cosmology and a set of concentration
     tables must be specified */
#ifdef TESTMODE
  printf("\n");

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"Omega0%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  Omega0 = atof(buffer2);
  printf("    Omega0                      %g\n", Omega0);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"OmegaLambda%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  OmegaLambda = atof(buffer2);
  printf("    OmegaLambda                 %g\n", OmegaLambda);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitMass_in_g%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  UnitMass_in_g = atof(buffer2);
  printf("    UnitMass_in_g               %g\n", UnitMass_in_g);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitLength_in_cm%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  UnitLength_in_cm = atof(buffer2);
  printf("    UnitLength_in_cm            %g\n", UnitLength_in_cm);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"UnitVelocity_in_cm_per_s%s",buffer2);
  checkforerror(error,elements,&buffer[0]);   
  UnitVelocity_in_cm_per_s = atof(buffer2);
  printf("    UnitVelocity_in_cm_per_s    %g\n", UnitLength_in_cm);

  error = fgets(buffer,sizeof(buffer),fp);
  elements = sscanf(buffer,"NameConcentrationTables%s",
		    NameConcentrationTables);
  checkforerror(error,elements,&buffer[0]);
  printf("    NameConcentrationTables     %s\n", NameConcentrationTables);

  UnitMass_in_Msun         = UnitMass_in_g/SOLARMASSG;
  UnitLength_in_kpc        = UnitLength_in_cm*1.0e3/MPCINCM;
  UnitVelocity_in_km_per_s = UnitVelocity_in_cm_per_s/1.0e5;
#endif

#ifndef TESTMODE
  if (RdTablesOn == 0 && VcTablesOn == 0) {
    fprintf(stderr, "\nError: no tables selected for output!\n");
    fprintf(stderr, "Check parameters file '%s' and try again - Exit\n",
	    filename); fflush(stderr);
    exit(EXIT_FAILURE);
  }
#endif
  printf("\n");

  return;
}
  
 
void checkforerror(char *error, int elements, char *buf)
{
  char	buffer2[256]; 
  
  if (error==NULL) {
    fprintf(stderr,"Error: couldn't read from parameter file. Stop.\n");
    exit(EXIT_FAILURE);
  }
  
  if (elements==0) {
    strncpy(buffer2,buf,strlen(buf));
    /* strncpy does not add a final \0 terminator. It has to be added by hand*/
    buffer2[strlen(buf)-1]='\0';
    fprintf(stderr,"Error: couldn't convert '%s'. Stop.\n",buffer2);
    exit(EXIT_FAILURE);
  }
}
