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
#include "LutefiskDefinitions.h"

/*	Global variables that have been declared in LutefiskGlobals.h */	

tmsms  msms;   /* Structure to hold data about the CID file. */

tionWeights gWeightedIonValues;

tParam gParam; /* Structure for the Lutefisk parameters loaded from Lutefisk.param */

char gSingAA[AMINO_ACID_NUMBER];

REAL_4 gMonoMass[AMINO_ACID_NUMBER];

REAL_4 gAvMass[AMINO_ACID_NUMBER];

INT_4 gNomMass[AMINO_ACID_NUMBER];

REAL_4 gWrongXCorrScore[WRONG_SEQ_NUM + 1];

REAL_4 gWrongIntScore[WRONG_SEQ_NUM + 1];

REAL_4 gWrongProbScore[WRONG_SEQ_NUM + 1];

REAL_4 gWrongQualityScore[WRONG_SEQ_NUM + 1];

REAL_4 gWrongComboScore[WRONG_SEQ_NUM + 1];

INT_4 gWrongIndex = 0;

INT_4 gTagLength;

INT_4 gSingleAACleavageSites = 0;

BOOLEAN gCorrectMass;

BOOLEAN gFirstTimeThru;

REAL_4 gElementMass[ELEMENT_NUMBER] = {
	 1.007825035, 		/* H */
	12.00,     		 	/* C */
	14.003074002, 	 	/* N */
	15.99491463, 	 	/* O */
	30.973762, 	 		/* P */
	31.972070698		/* S */
};

INT_4 gAminoAcidNumber = AMINO_ACID_NUMBER;

INT_4 gMonoMass_x100[AMINO_ACID_NUMBER];	/*Values assigned in 
											CreateGlobalIntegerMassArrays*/
INT_4 gElementMass_x100[ELEMENT_NUMBER];	/*Values assigned in 
											CreateGlobalIntegerMassArrays*/
INT_4 gMultiplier;							/*Value assigned in 
											CreateGlobalIntegerMassArrays*/
INT_4 gNodeCorrection[AMINO_ACID_NUMBER * AMINO_ACID_NUMBER];	/*Value assigned in 
											CreateGlobalIntegerMassArrays*/
INT_4 gElementCorrection[ELEMENT_NUMBER]; 	/*Value assigned in 
											CreateGlobalIntegerMassArrays*/
INT_4 gAvMonoTransition, gWater, gAmmonia, gCO, gAvResidueMass, gGraphLength;
											
INT_4 gGapList[MAX_GAPLIST];
INT_4 gGapListIndex;
INT_4 gEdmanData[MAX_PEPTIDE_LENGTH][AMINO_ACID_NUMBER];
INT_4 gMaxCycleNum;

REAL_4 H2O;
REAL_4 NH3;

REAL_4 gIonTypeWeightingTotal;

BOOLEAN	gDatabaseSeqCorrect;
