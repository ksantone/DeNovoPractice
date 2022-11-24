/*********************************************************************************************
Lutefisk is software for de novo sequencing of peptides from tandem mass spectra.
Copyright (C) 1995  Richard S. Johnson

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

Contact:

Richard S Johnson
4650 Forest Ave SE
Mercer Island, WA 98040

jsrichar@alum.mit.edu
*********************************************************************************************/

#ifndef _LUTEFISK_DEFS_
#define _LUTEFISK_DEFS_


#include <time.h>

#define DEBUG   /* Uncomment this line for debugging */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

#ifdef __LITTLE_ENDIAN 
#undef __LITTLE_ENDIAN
#endif
#ifdef __BIG_ENDIAN 
#undef __BIG_ENDIAN
#endif

#if (defined __MWERKS__ && __dest_os == __mac_os) || defined __OS_X
#define __BIG_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

#if TRUE//defined __MWERKS__ && __dest_os == __win32_os
#define __LITTLE_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

#if defined __SOLARIS
#define __BIG_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

#if defined __ALPHA
#define __LITTLE_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

#if defined __IRIX
#define __BIG_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

#if defined __AIX
#define __BIG_ENDIAN
#define CHAR    char
#define SCHAR   signed char
#define INT_2   signed short int
#define UINT_2  unsigned short int
#define INT_4   signed int
#define UINT_4  unsigned int
#define REAL_4  float
#define REAL_8  double
#define BOOLEAN unsigned char
#endif

#if defined __LINUX
#define __LITTLE_ENDIAN
#define CHAR	char
#define SCHAR	signed char
#define	INT_2	signed short int
#define UINT_2	unsigned short int
#define	INT_4	signed int
#define	UINT_4	unsigned int
#define	REAL_4	float
#define	REAL_8	double
#define BOOLEAN	unsigned char
#endif

       
/*	Define a few constants.	*/
#define AMINO_ACID_NUMBER 25	/*Number of amino acids.*/

/* AMINO ACIDS */
#define A  0
#define R  1
#define N  2
#define D  3
#define C  4
#define E  5
#define Q  6
#define G  7
#define H  8
#define I  9
#define L 10
#define K 11
#define M 12
#define F 13
#define P 14
#define S 15
#define T 16
#define W 17
#define Y 18
#define V 19


#define ELEMENT_NUMBER 6		/*Number of elements - H, C, N, O, P, and S.*/

#define HYDROGEN	0
#define CARBON		1
#define NITROGEN	2
#define OXYGEN		3
#define PHOSPHORUS  4
#define SULFUR		5


#define WATER ((gElementMass[HYDROGEN] * 2) + gElementMass[OXYGEN])	/*Mass of water.*/
#define AMMONIA (gElementMass[NITROGEN] + (gElementMass[HYDROGEN] * 3))	/*Mass of ammonia.*/
#define CO (gElementMass[CARBON] + gElementMass[OXYGEN])	/*Mass of carbon monoxide.*/
#define MAX_ION_NUM 500		/*500 Maximum number of fragment ions used in the final scoring.*/
#define MAX_DATA_POINTS_PER_GROUP 2500	/*The max number of data points per group of ions
										that exceed an ion threshold.  Each group is passed on
										to a function that sorts the ions.*/
#define SPECTRAL_WINDOW_WIDTH 120	/*Width of spectrum that has a limit on the number of ions.*/
#define MULTIPLIER_SWITCH 2.5	/*Used in determining gMultiplier*/
#define GRAPH_LENGTH 400000	/*The maximum graph size (ie, peptides must be less than 4000 Da).*/
#define AV_TO_MONO	0.999371395	/*The weighted average ratio between average and monoisotopic 
								amino acid masses.*/
#define MONO_TO_NOMINAL 0.99949725	/*The weighted average ratio between monoisotopic and
									nominal masses.*/
#define MONO_TO_AV 1.000629	/*To convert from monoisotopic to average mass.*/
#define NOMINAL_TO_MONO 1.000503003	/*To convert from nominal to monoisotopic mass.*/
#define NOMINAL_TO_AV 1.001132319	/*To convert from nominal to average mass.*/
#define C_NODE_VALUE 10	/*The value placed at the C-terminal nodes.*/
#define N_NODE_VALUE 10	/*The value placed at the N-terminal nodes.*/
#define ONE_EDGE_NODE_MAX  4000	/*Max number of one edged nodes.*/
#define MIN_MASS_PER_CHARGE 300	/*For charge states > 1 there is a minimum amount of mass 
								required to hold that charge.*/
#define MAX_PEPTIDE_LENGTH 60 	/*35 Maximum peptide length.*/
#define MAX_GAPLIST (AMINO_ACID_NUMBER * AMINO_ACID_NUMBER)
#define AV_RESIDUE_MASS 119	/*This is the weighted average amino acid residue mass.*/
#define TAG_MULTIPLIER 10
#define AV_MONO_TRANSITION	400	/*The factor that converts average to mono is transitioned over
								this mass range below the gParam.monoToAv value.*/
#define OVER_USED_IONS 0.9	/*For ions that are used for both y and b ions, the ionFound value
							is attenuated by this much.*/

/*The following affects the node values and/or subsequencing.*/
#define EDGE_EDGE_PENALTY 0.9		/*From a N-terminal one-edge node to a C-terminal
										one-edge node.*/
#define PROLINE_PENALTY 0.75	/*For gaps that might contain proline.*/
#define PRECURSOR_PENALTY 0.65	/*For gaps that encompass the precursor ion.*/
#define GLYCINE_PENALTY	0.5		/*For gaps that might contain glycine.*/
#define NODE_EDGE_PENALTY	0.4	/*From node to a C-terminal one-edge node.*/
#define NODE_NODE_PENALTY	0.2		/*These are multiplied against an extension score for an
									extension that uses two amino acids.*/
#define TOTALIONVAL_MULTIPLIER 1	/*Effects the node values that connect to C-terminus.*/
#define HIGH_CHARGE_MULT	0.5	/*Effects the node value if fragment has high charge (0-1).*/
#define HIGH_MASS_B_MULT	0.5	/*Effects the node value if b ions are more than precursor m/z (0-1).*/
#define HIGH_MASS_A_MULT	0.1	/*Effects the node value if a ions are more than 350 Da (0-1).*/

/*	These weighting values are used to weight the importance of various scoring parameters 
	when calculating the final score.  Used in LutefiskScore.*/
#define ATTENUATION_WEIGHT 0	/*1Presence of b and y ions.*/
#define INTENSITY_WEIGHT 1	/*2Ion current accounted for.*/
#define PEAKS_WEIGHT 0	/*0Peaks per residue divided by peaks per average residue.*/
#define NUMBER_WEIGHT 0	/*1Number of ions accounted for.*/
#define ATT_INT_PEAKS_NUM (ATTENUATION_WEIGHT + INTENSITY_WEIGHT + PEAKS_WEIGHT + NUMBER_WEIGHT)
#define INT_NUM_WEIGHT (NUMBER_WEIGHT + INTENSITY_WEIGHT)
#define INT_ATT_WEIGHT (INTENSITY_WEIGHT + ATTENUATION_WEIGHT)
#define INT_PEAKS_WEIGHT (INTENSITY_WEIGHT + PEAKS_WEIGHT)
#define MAX_X_CORR_NUM 5000	/*The maximum number of sequences to cross-corr score.*/
#define NEUTRAL_LOSS_MULTIPLIER 1	/*0.85 Fraction of ion signal for neutral loss ions (y-17, b-17,
									etc) that is counted in intensity score. Must be 0-1.*/
#define INTERNAL_FRAG_MULTIPLIER 1	/*0.5 Fraction of ion signal for internal fragment ions
										that is counted in intensity score. Must be 0-1.*/
#define HIGH_MASS_B_ION_MULTIPLIER 1	/*0.3 Fraction of ion signal for certain high mass b ions
										that is counted in intensity score. Must be 0-1.*/
#define HIGH_MASS_A_ION_MULTIPLIER 1	/*0.1Fraction of ion signal for certain high mass a ions
										that is counted in intensity score. Must be 0-1.*/
#define HIGH_CHARGE_Y_ION_MULTIPLIER 1	/*0.65 Fraction of ion signal for y ions that have the same
										charge as the precursor. Must be 0-1.*/
#define TWO_AA_EXTENSION_MULTIPLIER 1	/*1Fraction of ion signal for two amino acid extensions.
										Must be 0-1.*/
#define OXMET_MULTIPLIER 1				/*1Fraction of ion signal for losses of 46 u when mass
										accuracy is sufficient to differentiate oxMet.*/
#define PHE_MULTIPLIER 1				/*0.5Fraction of ion signal for losse of 46 u when mass
										accuracy is not sufficient to differentiate oxMet.*/
										
#define SIZEOF_SPECTRA_SMALL 2048	/* Arrays for FFT analysis must be a power of 2 in size */
#define SIZEOF_SPECTRA_BIG   4096

#define TAG_CUTOFF	50		/*Percentage of total ion current in tag region that must be
							accounted for before the sequence is considered a tag.*/
#define GOLDEN_BOY_CUTOFF	0.55	/*Fraction of average goldenBoy intensity as cutoff.*/
#define	GOLDEN_BOY_MAX		600		/*Golden boys cannot be greater than this m/z*/
#define IONFOUND_ISOTOPE	0.5	/*Value in ionFound for positions that could be isotopes of
								possible b or y ions.*/
#define SIGNAL_NOISE	4	/*This is the signal to noise required for an ion to be used.*/
#define MAX_QTOF_SEQUENCES	150	/*Max number of sequences for scoring using qtofErr*/
#define MAX_DATABASE_SEQ_NUM	50	/*Max number of sequences derived from database match*/
#define WRONG_SEQ_NUM	100	/*Number of wrong masses (and sequences) for comparing to correct*/



/* MACROS */
#define closeEnough(x,tol)       (((x)<0)? (-(x)<=tol):(x)<=tol)
#define MolWeightOf_Pos(mz,z)    ((mz - monoisotopicElementMass[HYDROGEN]) * (REAL_4)(z))
#define MolWeightOf_Neg(mz,z)    ((mz + monoisotopicElementMass[HYDROGEN]) * (REAL_4)(z))
#define mzOf_Pos(MW,z)           (((MW)/(REAL_4)(z)) + monoisotopicElementMass[HYDROGEN])
#define mzOf_Neg(MW,z)           (((MW)/(REAL_4)(z)) - monoisotopicElementMass[HYDROGEN])
#define chargeOf_Pos(MW,mz)      ((INT_4)(MW)/(mz -  monoisotopicElementMass[HYDROGEN]) + 0.1)   
#define SWAP(a,b)                tempr=(a);(a)=(b);(b)=tempr



/*	Define a few structs.	*/
extern struct MSData 		/*Used to hold ion information.*/
{
	REAL_4 mOverZ;
	INT_4 intensity;
	INT_4 normIntensity;
	struct MSData *next;
} MSData;

typedef struct {
	REAL_4 mOverZ;
	INT_4  intensity;
}tRawMSData;

typedef struct
{
	INT_4      numObjects;
	INT_4      limit;
	INT_4      sizeofobject;
	INT_4      growNum;
	tRawMSData   *mass;
}tRawMSDataList;

typedef struct {
	INT_4  index;
	REAL_4 mOverZ;
	INT_4 intensity;
	INT_4 normIntensity;
}tMSData;

typedef struct
{
	INT_4      numObjects;
	INT_4      limit;
	INT_4      sizeofobject;
	INT_4      growNum;
	tMSData   *mass;
}tMSDataList;

typedef struct 	/*Structure to hold data about the CID file.*/
{
	REAL_4 scanMassLow;
	REAL_4 scanMassHigh;
	char  instrument[10];
	char  centroidOrProfile;
}tmsms;

extern tmsms msms;

typedef struct
{
	INT_4 a;
	INT_4 a_minus17or18;
	INT_4 b;
	INT_4 b_minus17or18;
	INT_4 b_minusOH;
	INT_4 b_minusOH_minus17;
	INT_4 c;
	INT_4 d;
	INT_4 v;
	INT_4 w;
	INT_4 x;
	INT_4 y;
	INT_4 y_minus2;
	INT_4 y_minus17or18;
	INT_4 z_plus1;
}tionWeights;

extern tionWeights gWeightedIonValues;


typedef struct 				/*Used to hold Lutefisk's parameters.*/
{
	char		fMonitor;
	char		fVerbose;
	clock_t 	startTicks;
	clock_t 	searchTime;
	char		paramFile[256];
	char		outputFile[256];
	char 		cidFilename[256];
	char            detailsFilename[256];
	char            residuesFilename[256];
	REAL_4 		peptideMW;
	INT_4 		chargeState;
	BOOLEAN		maxent3;
	REAL_4 		fragmentErr;
	REAL_4		qtofErr;
	REAL_4 		ionOffset;
	char 		fragmentPattern;
	REAL_4 		cysMW;
	char 		proteolysis;
	char 		centroidOrProfile;
	INT_4 		monoToAv;
	REAL_4 		ionsPerWindow;
	char 		aaPresent[AMINO_ACID_NUMBER];
	char 		aaAbsent[AMINO_ACID_NUMBER];
	REAL_4 		modifiedNTerm; 
	REAL_4 		modifiedCTerm; 
	REAL_4 		tagNMass;
	INT_4		shoeSize;
	char 		tagSequence[MAX_PEPTIDE_LENGTH];
	REAL_4 		tagCMass; 
	REAL_4		outputThreshold;
	INT_4 		finalSeqNum;
	INT_4 		topSeqNum;
	REAL_4 		extThresh;
	INT_4 		maxExtNum;
	INT_4 		maxGapNum;
	INT_4		outputSeqNum;
	REAL_4 		peakWidth;
	REAL_4 		ionThreshold; /* Fraction */
	INT_4       intThreshold; /* Actual intensity threshold */
	REAL_4 		peptideErr;
	BOOLEAN	 	autoTag;
	BOOLEAN 	edmanPresent;
	char		edmanFilename[256];
	REAL_4 		ionsPerResidue;
	char		CIDfileType;
	char		databaseSequences[256];
	BOOLEAN		quality;
	INT_4		wrongSeqNum;
	
	INT_4		topSeqNum_orig;
	REAL_4 		peptideMW_orig;
	REAL_4 		peptideErr_orig;
	REAL_4 		fragmentErr_orig;
	REAL_4 		peakWidth_orig;
	INT_4 		monoToAv_orig;
	REAL_4		qtofErr_orig;
	REAL_4 		ionOffset_orig;
	REAL_4 		cysMW_orig;
	REAL_4 		tagNMass_orig;
	REAL_4 		tagCMass_orig; 
	REAL_4		modifiedNTerm_orig;
	REAL_4		modifiedCTerm_orig;
        INT_4 		maxGapNum_orig;
	char		outputFile_orig[256];
        
}tParam;

extern tParam gParam;


extern struct Sequence		/*Used to hold sequence info during subsequencing.*/
{
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	INT_4 peptideLength;
	INT_4 score;
	INT_4 nodeValue;
	INT_2 nodeCorrection;
	INT_4 gapNum;
	struct Sequence *next;
} Sequence;


typedef struct
{
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	INT_4 peptideLength;
	INT_4  score;
	INT_4 nodeValue;
	INT_4 gapNum;
}tsequence;

typedef struct
{
	INT_4      numObjects;
	INT_4      limit;
	INT_4      sizeofobject;
	INT_4      growNum;
	tsequence *seq;
}tsequenceList;

extern struct SequenceScore		/*Used to hold sequence info during final scoring procedure.*/
{
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	REAL_4 intensityScore;
	REAL_4 intensityOnlyScore;
	REAL_4 crossDressingScore;
	REAL_4 stDevErr;
	REAL_4 calFactor;
	REAL_4 quality;	/*new 6/20/01 this is contiguous single aa divided by actual peptide length*/
	REAL_4 length;	/*new 6/20/01 this is the unfudged unadulterated peptide length, counting
					dipeptide residues as two amino acids*/
	REAL_8 probScore;	/*7/30/03 added Pevzner probability score*/
	REAL_4 comboScore;	/*final combined score*/
	INT_4 cleavageSites;
	INT_2 peptideSequence[MAX_PEPTIDE_LENGTH];
	char databaseSeq;
	INT_4 rank;
	struct SequenceScore *next;
} SequenceScore;

typedef struct extension {
	char	gapSize;
	BOOLEAN	singleAAFLAG;
	INT_4	mass;
	INT_4	score;
	INT_2   nodeCorrection;
} extension;


extern char gSingAA[AMINO_ACID_NUMBER];
  
extern REAL_4 gMonoMass[AMINO_ACID_NUMBER];

extern REAL_4 gAvMass[AMINO_ACID_NUMBER];

extern INT_4 gNomMass[AMINO_ACID_NUMBER];

extern REAL_4 gElementMass[ELEMENT_NUMBER];

extern INT_4 gGapList[MAX_GAPLIST];

extern INT_4 gAminoAcidNumber;

extern INT_4 gGapListIndex;

extern INT_4 gEdmanData[MAX_PEPTIDE_LENGTH][AMINO_ACID_NUMBER], gMaxCycleNum;

extern REAL_4 H2O, NH3;

extern INT_4 gElementMass_x100[ELEMENT_NUMBER];	/*Values assigned in 
												CreateGlobalIntegerMassArrays*/
extern INT_4 gMonoMass_x100[AMINO_ACID_NUMBER];	/*Values assigned in 
												CreateGlobalIntegerMassArrays*/
extern INT_4 gMultiplier;	/*Value determined in CreateGlobalIntegerMassArrays*/

extern INT_4 gNodeCorrection[MAX_GAPLIST]; /*Value determined in CreateGlobalIntegerMassArrays*/

extern INT_4 gElementCorrection[ELEMENT_NUMBER]; /*Value determined in CreateGlobalIntegerMassArrays*/

extern INT_4 gAvMonoTransition, gWater, gAmmonia, gCO, gAvResidueMass, gGraphLength;

extern REAL_4 gWrongXCorrScore[WRONG_SEQ_NUM + 1];

extern REAL_4 gWrongIntScore[WRONG_SEQ_NUM + 1];

extern REAL_4 gWrongProbScore[WRONG_SEQ_NUM + 1];

extern REAL_4 gWrongQualityScore[WRONG_SEQ_NUM + 1];

extern REAL_4 gWrongComboScore[WRONG_SEQ_NUM + 1];

extern INT_4 gSingleAACleavageSites;

extern INT_4 gWrongIndex;

extern INT_4 gTagLength;

extern BOOLEAN gCorrectMass;

extern BOOLEAN gFirstTimeThru;

extern REAL_4 gIonTypeWeightingTotal;

extern BOOLEAN gDatabaseSeqCorrect;

#endif /* _LUTEFISK_DEFS_ */
