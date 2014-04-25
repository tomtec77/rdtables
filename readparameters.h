/**
 * @file readparameters.h
 * @brief Declarations of the global variables that store the parameters
 * given as input to the program
 * @author Tomas E. Tecce
 */
extern char Path1[256];

extern double BaryonFrac;

#ifdef TESTMODE
extern double Omega0, OmegaLambda;
extern double UnitMass_in_g, UnitLength_in_cm, UnitVelocity_in_cm_per_s;
extern double UnitMass_in_Msun, UnitVelocity_in_km_per_s, UnitLength_in_kpc;
#endif

extern int NumOutputTimes;

extern double TimeOfFinalOutput;
extern double TimeBetOutputs;

extern int NumCbins;

extern int NumOutputTimes;

extern double TimeOfFinalOutput;
extern double TimeBetOutputs;

extern int NumMdiscBins;

extern double StartingMdisc;
extern double WidthOfMdiscBins;

extern int NumLambdaBins;

extern double StartingLambda;
extern double WidthOfLambdaBins;

extern int NumMbulgeBins;

extern double StartingMbulge;
extern double WidthOfMbulgeBins;

/*extern double MassFractionForVc;*/
extern double OpticalRadius;

extern unsigned int OutputListOn;
extern unsigned int ConcentrationListOn;
extern unsigned int MbulgeListOn;
extern unsigned int BulgeAngularMomentumOn;
extern unsigned int RdTablesOn, VcTablesOn;

extern char NameOutputList[256];
extern char NameMbulgeList[256];
extern char NameConcentrationTables[256];
