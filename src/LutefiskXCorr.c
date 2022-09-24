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

/*
	Richard S. Johnson
	6/96 - ?
	
		LutefiskXP is a program designed to aid in the interpretation of CID data of peptides.  
	The main assumptions are that the data is of reasonable quality, the N- and C-terminal
	modifications (if any) are known, and the precursor ion charge (and therefore the 
	peptide molecular weight) are known.  The ultimate goal here is to develop code that
	can utilize msms data in conjunction with ambiguous and incomplete Edman sequencing data,
	sequence tags, peptide derivatization, and protein or est database searches.  An older 
	version of LutefiskXP has been written in FORTRAN and runs on 68K Macs that have an fpu
	(1991, 39th ASMS Conference on Mass Spectrometry and Allied Topics, Nashville, TN, pp 1233-
	1234).  This is a different and improved algorithm partly inspired by Fernandez-de-Cossjo, 
	et al. (1995) CABIOS Vol. 11 No. 4 pp 427-434.  Combining this msms interpretation algorithm
	with Edman sequencing, database searches, and derivatization is entirely of my own design;
	J. Alex Taylor implemented the changes in the FASTA code (Bill Pearson, U. of VA) so that
	the LutefiskXP output can be read directly by the modified FASTA program.  In addition, there
	were a number of additional critical changes made to FASTA to make it more compatible with
	msms sequencing data.
	
		The trademark LutefiskXP was chosen at random, and is not meant to imply any similarity
	between this computer program and the partially base-hydrolzyed cod fish of the same name
	(minus XP). 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Lutefisk headers*/
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"


#define MAXPARENTCHARGE 9
#define SIDE_PEAK_ATT	0.75	/*Peak heights on sides of main peak in mock spectrum.*/
#define PLUS1_NEUT_LOSS_ATT	0.5	/*Peak height for neutral losses of singly charged precursors.*/
#define NEUT_LOSS_ATT 0.1	/*Peak height for neutral losses of multiply charged precursors.*/
#define BAD_B_ATT	0.05	/*Peak heights for b ions that are not very likely.*/
#define A_ATT	0.5	/*Peak heights for a ions.*/
#define BAD_A_ATT 0.05	/*Peak heights for a ions that are not very likely.*/
#define INT_FRAG_ATT 0.1	/*peak heights for internal fragment ions.*/
#define BAD_Y_ATT	0.05	/*Peak heights for y ions that are not very likely.*/

extern REAL_4 	*spectrum1;
extern REAL_4 	*spectrum2;
extern REAL_4 	*tau;

UINT_4 SIZEOF_SPECTRA;	/* Smallest power of 2 (for cross-correlation)*/
REAL_4 gSidePeakAtt = SIDE_PEAK_ATT;


/*Here's some globals that are specific to this file.  They are two amino acid nominal masses
*times 100 for Cys, Arg, His, and Lys.  These get modified at the start of the ScoreSequences
*function in order to accomodate different alkyl groups on cysteine.
*/
INT_4 gCysPlusXCorr[AMINO_ACID_NUMBER] = {
	174, 259, 217, 218, 206, 232, 231, 160, 240, 216,
	216, 231, 234, 250, 200, 190, 204, 289, 266, 202,
	0,0,0,0,0
};

INT_4 gArgPlusXCorr[AMINO_ACID_NUMBER] = {
	227, 312, 270, 271, 259, 285, 284, 213, 293, 269,
	269, 284, 287, 303, 253, 243, 257, 342, 319, 255,
	0,0,0,0,0
};

INT_4 gHisPlusXCorr[AMINO_ACID_NUMBER] = {
	208, 293, 251, 252, 240, 266, 265, 194, 274, 250,
	250, 265, 268, 284, 234, 224, 238, 323, 300, 236,
	0,0,0,0,0
};

INT_4 gLysPlusXCorr[AMINO_ACID_NUMBER] = {
	199, 284, 242, 243, 231, 257, 256, 185, 265, 241,
	241, 256, 259, 275, 225, 215, 229, 314, 291, 227,
	0,0,0,0,0
};

INT_4 lowMassIonMass[AMINO_ACID_NUMBER] = {
	44, 112, 87, 88, 76, 102, 102, 30, 110, 86, 
	86, 129, 104, 120, 70, 60, 74, 159, 136, 72,
	0,0,0,0,0
};

REAL_4 lowMassIonIntFactor[AMINO_ACID_NUMBER] = {
	0.1, 0.1, 0.1, 0.1, 0.0, 0.1, 0.1, 0.0, 0.5, 0.5, 
	0.5, 0.25, 0.1, 0.5, 0.3, 0.1, 0.1, 0.2, 0.3, 0.2,
	0,0,0,0,0
};

extern void CrossCorrScoreTheSeq(struct SequenceScore *currScorePtr);

/***************************YXCorrCalc*****************************************************
*
*	This function calculates the mass of singly charged y ions for a given cleavage,
*	and a given sequence.  It returns a REAL_4.
*
*/
REAL_4 YXCorrCalc(INT_4 i, struct SequenceScore *currScorePtr, REAL_4 YionStart, 
					INT_4 seqLength)
{
	REAL_4 Yion;
	INT_4 j, k;
	char test;
	REAL_8 mToAFactor;
	
	Yion = YionStart;
	for(j = i; j < seqLength; j++)
	{
		test = TRUE;
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(currScorePtr->peptide[j] == gNomMass[k])
			{
				test = FALSE;
				Yion += gMonoMass[k];
				break;
			}
		}
		if(test)	/*its a 2 aa extension*/
		{
			Yion = Yion + ((currScorePtr->peptide[j])) + 0.07;	/*0.07 is a guess at the mass defect*/
		}
	}
	
	if(Yion > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(Yion >= (gParam.monoToAv - AV_MONO_TRANSITION))
		{
			mToAFactor = (gParam.monoToAv - Yion) / AV_MONO_TRANSITION;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(Yion >= (gParam.monoToAv - AV_MONO_TRANSITION))
	{
		Yion = Yion * mToAFactor;
	}


	return(Yion);
}
/******************************BXCorrCalc***************************************************
*
*	This function calculates the mass of singly charged b ions for a given cleavage,
*	and a given sequence.  It returns a REAL_4.
*
*/
REAL_4 BXCorrCalc(INT_4 i, struct SequenceScore *currScorePtr, REAL_4 BionStart)
{
	REAL_4 Bion;
	INT_4 j, k;
	char test;
	REAL_8 mToAFactor;
	
	Bion = BionStart;
	for(j = 0; j < i; j++)
	{
		test = TRUE;
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(currScorePtr->peptide[j] == gNomMass[k])
			{
				test = FALSE;
				Bion += gMonoMass[k];
				break;
			}
		}
		if(test)	/*its a 2 aa extension*/
		{
			Bion = Bion + ((currScorePtr->peptide[j])) + 0.07;	/*0.07 is a guess at the mass defect*/
		}
	}
	
	if(Bion > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(Bion >= (gParam.monoToAv - AV_MONO_TRANSITION))
		{
			mToAFactor = (gParam.monoToAv - Bion) / AV_MONO_TRANSITION;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(Bion >= (gParam.monoToAv - AV_MONO_TRANSITION))
	{
		Bion = Bion * mToAFactor;
	}
		
	return(Bion);
}


/*********************************FindNChargeXCorr********************************************
*
*	Counts the number of charged residues in the sequence.
*
*/
INT_4 FindNChargeXCorr(struct SequenceScore *currScorePtr)
{
	INT_4 j, k, seqLength, nChargeCount, *pSeq;
	
/*	Initialize.*/
	nChargeCount = 1;
	seqLength = 0;

/*
	Determine the sequence length.
*/
	pSeq = &currScorePtr->peptide[0];
	while(*pSeq != NULL)
	{
		seqLength++;
		pSeq++;
	}
/*
	Find one aa extensions that are either R, H, or K.
*/
	for(j = 0; j < seqLength; j++)
	{
		if((currScorePtr->peptide[j] == gNomMass[R]) || (currScorePtr->peptide[j] == gNomMass[H]) || 
			(currScorePtr->peptide[j] == gNomMass[K]))
		{
			nChargeCount += 1;
		}
	}
/*
	Here I look for two amino acid extensions containing Arg, His, or Lys.
*/
	for(j = 0; j < seqLength; j++)	
	{
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(currScorePtr->peptide[j] == gArgPlusXCorr[k] || 
				currScorePtr->peptide[j] == gHisPlusXCorr[k] || 
				currScorePtr->peptide[j] == gLysPlusXCorr[k])
			{
				nChargeCount += 1;
			}
		}
	}
	return(nChargeCount);
}


/**************************	DoCrossCorrelationScoring*******************************************
*
*
*/
void DoCrossCorrelationScoring(struct SequenceScore *firstScorePtr, struct MSData *firstMassPtr) 
{

	INT_4 		i, seqNum;						
	REAL_4		normalizedScore;		/* The max absolute cross-correlation value
										*   (Later used as the normalizing factor)*/
	REAL_4		tauDiff, intensityAccountedFor, autocorrelation;
	struct SequenceScore *currSeqPtr;


	if (!firstScorePtr) 
	{
		return;
	}
	
	/*For high accuracy data, make peaks w/ no extra width*/
	if(gParam.qtofErr != 0)
	{
		if(gParam.qtofErr < 0.25)
		{
			gSidePeakAtt = 0;
		}
		else
		{
			gSidePeakAtt = SIDE_PEAK_ATT;
		}
	}
	else
	{	
		gSidePeakAtt = SIDE_PEAK_ATT;
	}

/*Do the funny normalization of intensities that Eng does in Sequest.*/
	CalcNormalizedExptPeaks(firstMassPtr);

/*Set aside memory for the spectra that are cross-correlated.*/
	SetupCrossCorrelation();

	if (!spectrum1 || !spectrum2 || !tau) 
	{
		/* Do not proceed with cross-correlation because we couldn't get the memory*/
		return;
	}
	else  
	{
		FillInSpectrum1(firstMassPtr);	/*Using the list of ions from firstMassPtr, generate
										a "real" spectrum for cross-correlating.  These mass
										values will be 2 times the actual number in order to
										add to the specificity of the cross-correlation.*/
	}
	
	/*Do an autocorrelation of the spectrum.  This is the normalizing factor used later*/
	
	CrossCorrelate(spectrum1-1, spectrum1-1, (UINT_4) SIZEOF_SPECTRA, tau-1);
	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		if(tau[i] < 1)
		{
			tau[i] = 0;
		}
	}
	
	
	intensityAccountedFor = 0.0;
	for(i = 1; i < 250; i++)
	{
		tauDiff = tau[i] - tau[SIZEOF_SPECTRA - i];
		if(tauDiff < 0)
		{
			tauDiff = tauDiff * -1;
		}
		intensityAccountedFor += tauDiff;
	}
	autocorrelation = tau[0] - intensityAccountedFor/250; /*should be equal to tau(0)*/

	
	/*debug spectrum1*/
	/*printf("data spectrum \n");	
	for(i = 0; i < SIZEOF_SPECTRA; i++)	
	{
		if(spectrum1[i] != 0)
		{
			printf("%f  %f \n",(REAL_4)i/2,spectrum1[i]);
		}
	}*/

/*Count the number of sequences.*/
	seqNum = 0;
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		seqNum++;
		currSeqPtr = currSeqPtr->next;
	}
	
/*If there are too many, then cut the number of sequences to be cross-correlated.*/
	if(seqNum > MAX_X_CORR_NUM)
	{
		seqNum = MAX_X_CORR_NUM;
	}
	
/*Cross-correlate the sequences.  Only do the top intensity-scorers.*/
	for(i = 1; i <= seqNum; i++)
	{	
		currSeqPtr = firstScorePtr;
		while(currSeqPtr != NULL)
		{
			if(i == currSeqPtr->rank)
			{
				CrossCorrScoreTheSeq(currSeqPtr);
			}
			currSeqPtr = currSeqPtr->next;
		}

	}
		
	/* Normalize the cross-correlation results to 1.0*/
	normalizedScore = 0.0 ;	/* First find the highest score*/
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		if (currSeqPtr->crossDressingScore > normalizedScore) 
		{
			normalizedScore = currSeqPtr->crossDressingScore;
		}
		currSeqPtr = currSeqPtr->next;	/* Point to the next struct in the linked list to 
										* continue the while loop.	*/
	}	
	/* Use the highest score to normalize them all*/
	if (normalizedScore > 0.0) 
	{
		normalizedScore = 1/normalizedScore;
	}
	
	/*normalizedScore = 1;*/	/*to keep from normalizing*/
	/*normalizedScore = 1 / gParam.peptideMW;*/ 
	if(autocorrelation == 0)
	{	
		printf("Avoiding divide by zero.");
                autocorrelation += 0.0001;  // JAT 2006.08.28
//		exit(1);
	}
	normalizedScore = 1 / autocorrelation;	/*another way to normalize*/
	
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		currSeqPtr->crossDressingScore = normalizedScore * currSeqPtr->crossDressingScore;

		currSeqPtr = currSeqPtr->next;	/* Point to the next struct in the linked list to 
										*  continue the while loop.	*/
	}	
	
	if (spectrum1)	free(spectrum1);	spectrum1 = NULL;
	if (tau)        free(tau);	          tau = NULL;
	
	return;
}
/**************************	CrossCorrScoreTheSeq*******************************************
*
*  This function calculates the cross-correlation score for ea. peptide passed to it.
*/
void CrossCorrScoreTheSeq(struct SequenceScore *currScorePtr) 
{

	INT_4 			i, j, k;					/* Loop indicies. */
	INT_4			chargeLimit;				/* Max # of daughter charges to consider */
	REAL_4 			intensityAccountedFor = 0.0;/* Sum of the intensities of the ions matched with fragment ions. */
	REAL_4 			Bion;						/* Mass of the singly charged B ion. */         
	REAL_4 			Yion;						/* Mass of the singly charged Y ion. */
	REAL_4			BionStart;					/* Mass of N-terminal group (H, Ac, etc)*/
	REAL_4 			YionStart;					/* Mass of Y-terminal group (free or amidated)*/
	REAL_4 			BionMass;					/* Mass of the B ion in a particular charge state. */
	REAL_4 			AionMass;					/* Mass of the A ion (B ion - CO) in a particular charge state. */
	REAL_4 			YionMass;					/* Mass of the Y ion in a particular charge state. */
	REAL_4			proAttenuation;				/* Attenuations intensity for y ions next to Pro*/
	REAL_4			glyAttenuation;				/* Attenuates intensity for y ions next to Gly*/
	REAL_4 			fragmentMass;				/* Mass of the internal fragment ion. */
	REAL_4 			parent;				 		/* The m/z of the parent ion. */
	REAL_4 			tolerance;					/* The daughter ion error tolerance. */
	REAL_4 			offset; 					/* The mass offset for daughter ions. */
	REAL_4 	parentMinH2O, parentMinNH3, parentMin2H2O, parentMin2NH3;
	REAL_4	tauDiff;
	REAL_4	highMassRange, lowMassRange;

	REAL_4 			fullIntensity;				/* Intensity used by cross-correlation when creating dummy spectra. */
	INT_4			seqLength, nChargeCount, cChargeCount;
	INT_4			*pSeq;
	char			NTermTest, widePeak;

	parent    = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	tolerance = gParam.fragmentErr;
	if(gParam.fragmentPattern == 'L')	/*ion masses outside of this range are not penalized if not found*/
	{
		lowMassRange = parent * 0.333;	/*so-called 1/3 rule*/
		highMassRange = 2000;	/*mass limit for Deca*/
	}
	else
	{
		lowMassRange = 146;	/*y1 for Lys*/
		highMassRange = 2 * parent;	/*often the very high mass ions are missing*/
	}
	offset    = 0;	/*the offset has already been applied*/
	chargeLimit = gParam.chargeState;
	if(gParam.fragmentErr <= 0.75)
	{
		widePeak = FALSE;	/*mock peak widths are 1.5 dalton*/
	}
	else
	{
		widePeak = TRUE;	/*mock peak widths are 2.5 daltons*/
	}
		
	if (spectrum2) {
		free(spectrum2);   /* Throw away the old data (if any exists) */
		spectrum2 = NULL;
	}	
	spectrum2 = (REAL_4 *) malloc(SIZEOF_SPECTRA*sizeof(REAL_4));
	if(spectrum2 == NULL)
	{
		printf("Out of memory");
		exit(1);
	}
	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		spectrum2[i] = 0;
	}
	if (!spectrum2) return;

	seqLength = 0;
	
/* 
	Determine the mass of the N and C-termini.
*/
	BionStart = gParam.modifiedNTerm;
	
	YionStart = gParam.modifiedCTerm + (2 * gElementMass[HYDROGEN]);

	
/*	
	Setup the values for nChargeCount and cChargeCount; FindNChargeXCorr also 
	removes the factor of 100 from the sequence masses.
*/
	cChargeCount = 1;
	nChargeCount = FindNChargeXCorr(currScorePtr);
	
/*
	Test if N-terminal amino acid is one or two residues.
*/
		NTermTest = TRUE;	
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(currScorePtr->peptide[0] == gNomMass[k])
			{	
				NTermTest = FALSE;	/*TRUE if a two amino acid step.*/
			}
		}
/*	
	Find the sequence length, so I can loop through it.
*/
	pSeq = &currScorePtr->peptide[0];
	seqLength = 0;
	while(*pSeq != NULL)
	{
		seqLength++;
		pSeq++;
	}

/*
	Here's where I loop through the sequence.
	=================================================================
*/
	for (i = seqLength - 1; i > 0; i--) 
	{
		fullIntensity = 50;	/*Arbitrary, but it should match the spectrum1 max value.*/

/*	Figure out what Bion and Yion should be.*/
		Bion = BXCorrCalc(i, currScorePtr, BionStart);
		Yion = YXCorrCalc(i, currScorePtr, YionStart, seqLength);

/*	Figure out what cChargeCount and nChargeCount should be.*/

		if((currScorePtr->peptide[i] == gNomMass[R]) || 
			(currScorePtr->peptide[i] == gNomMass[H]) || 
			(currScorePtr->peptide[i] == gNomMass[K]))
		{
			cChargeCount += 1;
			nChargeCount -= 1;
		}
		else	/*Check to see if its a two amino acid combo that could contain Arg, His, or Lys.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(currScorePtr->peptide[i] == gArgPlusXCorr[j] || 
					currScorePtr->peptide[i] == gHisPlusXCorr[j] || 
					currScorePtr->peptide[i] == gLysPlusXCorr[j])
				{
					cChargeCount += 1;
					nChargeCount -= 1;
					break;
				}
			}
		}
		
/* 
	Look for ea. possible charge state.>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/
		for (j = 1; j <= chargeLimit; j++) {   
			BionMass = (Bion + (j - 1))/j;
			AionMass = (Bion - CO + (j - 1))/j;
			YionMass = (Yion + (j - 1))/j;
		 	/* Once we are looking for fragment ions of charge 3 or above these should be considered
		 	*  less likely and so, if cross-correlation is on, these peaks will be given an intensity
		 	*  quarter of what it would be for a charge 1 or 2 fragment ion. */
		 	if (j == 3) fullIntensity = fullIntensity / 4;

		 	
		 	/* Dummy up a theoretical spectrum for the peptide.
			*  Maker sure there are chargeable amino acids for the given charge - j - and don't
			*  bother with b and a ions that have just a proton or have just a proton and one aa.*/
				
			if(nChargeCount >= j && i != seqLength && (i != 1 || NTermTest))
			{
				if((BionMass * j) > (j-1) * 500)	/*Make sure there is enough mass for the charge.*/
				{
	 				if (BionMass < msms.scanMassHigh && BionMass < highMassRange
	 					&& BionMass > lowMassRange) 
					{
						if(gParam.chargeState == 1)	/*If the precursor is singly charged.*/
						{
							BionMass = BionMass * 2;
							NH3 = NH3 * 2;
							H2O = H2O * 2;
							
							AddPeakToSpectrum(spectrum2, BionMass, fullIntensity);
							AddPeakToSpectrum(spectrum2, BionMass - NH3/j, fullIntensity * PLUS1_NEUT_LOSS_ATT);
							AddPeakToSpectrum(spectrum2, BionMass - H2O/j, fullIntensity * PLUS1_NEUT_LOSS_ATT);

							BionMass = BionMass * 0.5;
							NH3 = NH3 * 0.5;
							H2O = H2O * 0.5;
						}
						else	/*For multiply charged precursors.*/
						{
							/*If the +1 b ion has an m/z less than the precursor ion or if the number of chargeable
							amino acids in the +1 b ion equals the charge state of the precursor, then
							give these guys full intensity.*/
							if((j == 1 || j <= (gParam.chargeState - 1)) && 
								(BionMass < parent || nChargeCount == gParam.chargeState
								|| gParam.fragmentPattern == 'L'))
							{
								BionMass = BionMass * 2;
								NH3 = NH3 * 2;
								H2O = H2O * 2;
								
								AddPeakToSpectrum(spectrum2, BionMass, fullIntensity);
								AddPeakToSpectrum(spectrum2, BionMass - NH3/j, fullIntensity * NEUT_LOSS_ATT);
								AddPeakToSpectrum(spectrum2, BionMass - H2O/j, fullIntensity * NEUT_LOSS_ATT);
								
								BionMass = BionMass * 0.5;
								NH3 = NH3 * 0.5;
								H2O = H2O * 0.5;
							}
							else	/*Otherwise... (ie, >+1 b ions or +1 bions > precursor)*/
							{
								BionMass = BionMass * 2;
								NH3 = NH3 * 2;
								H2O = H2O * 2;
								
								AddPeakToSpectrum(spectrum2, BionMass, fullIntensity * BAD_B_ATT);
								AddPeakToSpectrum(spectrum2, BionMass - NH3/j, fullIntensity * BAD_B_ATT * NEUT_LOSS_ATT);
								AddPeakToSpectrum(spectrum2, BionMass - H2O/j, fullIntensity * BAD_B_ATT * NEUT_LOSS_ATT);
								
								BionMass = BionMass * 0.5;
								NH3 = NH3 * 0.5;
								H2O = H2O * 0.5;
							}
						}
	 				}
		
			 		if (AionMass < msms.scanMassHigh && AionMass < highMassRange && AionMass > lowMassRange) 
			 		{
						if(gParam.chargeState == 1)
						{
							AionMass = AionMass * 2;
							
							AddPeakToSpectrum(spectrum2, AionMass, fullIntensity * A_ATT);
							
			 				AionMass = AionMass * 0.5;
						}
						else
						{
/*	Give higher intensity for singly charged a2 ions.  Skip if LCQ data (?).*/
							if(j == 1 && ((NTermTest == TRUE && i == 1) || 
											(NTermTest == FALSE && i == 2)))
							{
								AionMass = AionMass * 2;
								
								AddPeakToSpectrum(spectrum2, AionMass, fullIntensity * A_ATT);
								
								AionMass = AionMass * 0.5;
							}
							else
							{
								if((j == 1 || j <= (gParam.chargeState - 1)) && 
									(BionMass < parent ||  nChargeCount == gParam.chargeState))
								{
									AionMass = AionMass * 2;
									
									AddPeakToSpectrum(spectrum2, AionMass, fullIntensity * A_ATT * BAD_A_ATT);
									
			 						AionMass = AionMass * 0.5;
								}
								else
								{
									AionMass = AionMass * 2;
									
									AddPeakToSpectrum(spectrum2, AionMass, fullIntensity * A_ATT * BAD_A_ATT * BAD_A_ATT);

			 						AionMass = AionMass * 0.5;
								}
							}
						}
			 		}	
				}
			}

			if(cChargeCount >= j && i != 0)
			{
				if((YionMass * j) > (j-1) * 500)
				{
	 				if (YionMass < msms.scanMassHigh && YionMass < highMassRange && YionMass > lowMassRange) 
					{
						proAttenuation = 1.0;
						glyAttenuation = 1.0;
						if(i > 2 && i < seqLength - 2)	/*don't attenuate at the ends*/
						{
							if(currScorePtr->peptide[i-1] == gNomMass[P])
							{
								proAttenuation = 0.2;	/*Attenuates y ions followin P*/
							}
							if(currScorePtr->peptide[i-1] == gNomMass[G])
							{
								/*glyAttenuation = 0.2;*/	/*Attenuates y ions followin G*/
							}
						}
						
						fullIntensity = fullIntensity * proAttenuation * glyAttenuation;
										
						if(gParam.chargeState == 1)
						{
							YionMass = YionMass * 2;
							NH3 = NH3 * 2;
							H2O = H2O * 2;
							
							AddPeakToSpectrum(spectrum2, YionMass, fullIntensity);
							AddPeakToSpectrum(spectrum2, YionMass - NH3/j, fullIntensity * PLUS1_NEUT_LOSS_ATT);
							AddPeakToSpectrum(spectrum2, YionMass - H2O/j, fullIntensity * PLUS1_NEUT_LOSS_ATT);

	 						YionMass = YionMass * 0.5;
	 						NH3 = NH3 * 0.5;
							H2O = H2O * 0.5;
						}
						else
						{	/*LCQ data has a preponderance of multiply charged y ions*/
							if(j == 1 || j <= (gParam.chargeState - 1) || 
								(i <= 2 && !NTermTest) || (i <= 1 && NTermTest) ||
								(i <= (INT_4)seqLength / 4 && gParam.fragmentPattern == 'L'))
							{
								YionMass = YionMass * 2;
	 							NH3 = NH3 * 2;
								H2O = H2O * 2;
								if(j == 1 || j <= (gParam.chargeState - 1) || 
								(i <= 2 && !NTermTest) || (i <= 1 && NTermTest))
								{
									AddPeakToSpectrum(spectrum2, YionMass, fullIntensity);	/*add some width to multiply charged*/
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - 1, fullIntensity);
										AddPeakToSpectrum(spectrum2, YionMass + 1, fullIntensity);
									}
									AddPeakToSpectrum(spectrum2, YionMass - NH3/j, fullIntensity * NEUT_LOSS_ATT);
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - NH3/j - 1, fullIntensity * NEUT_LOSS_ATT);
										AddPeakToSpectrum(spectrum2, YionMass - NH3/j + 1, fullIntensity * NEUT_LOSS_ATT);
									}
									
									AddPeakToSpectrum(spectrum2, YionMass - H2O/j, fullIntensity * NEUT_LOSS_ATT);
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - H2O/j - 1, fullIntensity * NEUT_LOSS_ATT);
										AddPeakToSpectrum(spectrum2, YionMass - H2O/j + 1, fullIntensity * NEUT_LOSS_ATT);
									}
								}
								else
								{
									AddPeakToSpectrum(spectrum2, YionMass, fullIntensity * 1/i);	/*add some width to multiply charged*/
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - 1, fullIntensity * 1/(i+ 1));
										AddPeakToSpectrum(spectrum2, YionMass + 1, fullIntensity * 1/(i+ 1));
									}
									
									AddPeakToSpectrum(spectrum2, YionMass - NH3/j, fullIntensity * NEUT_LOSS_ATT);
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - NH3/j - 1, fullIntensity * 1/(i+ 1) * NEUT_LOSS_ATT);
										AddPeakToSpectrum(spectrum2, YionMass - NH3/j + 1, fullIntensity * 1/(i+ 1) * NEUT_LOSS_ATT);
									}
									
									AddPeakToSpectrum(spectrum2, YionMass - H2O/j, fullIntensity * NEUT_LOSS_ATT);
									if(j > 1)
									{
										AddPeakToSpectrum(spectrum2, YionMass - H2O/j - 1, fullIntensity * 1/(i+ 1) * NEUT_LOSS_ATT);
										AddPeakToSpectrum(spectrum2, YionMass - H2O/j + 1, fullIntensity * 1/(i+ 1) * NEUT_LOSS_ATT);
									}
								}
								
								YionMass = YionMass * 0.5;
	 							NH3 = NH3 * 0.5;
								H2O = H2O * 0.5;
							}
							else
							{
								YionMass = YionMass * 2;
	 							NH3 = NH3 * 2;
								H2O = H2O * 2;
								
								AddPeakToSpectrum(spectrum2, YionMass, fullIntensity * BAD_Y_ATT);
								AddPeakToSpectrum(spectrum2, YionMass - NH3/j, fullIntensity * BAD_Y_ATT * NEUT_LOSS_ATT);
								AddPeakToSpectrum(spectrum2, YionMass - H2O/j, fullIntensity * BAD_Y_ATT * NEUT_LOSS_ATT);

								YionMass = YionMass * 0.5;
	 							NH3 = NH3 * 0.5;
								H2O = H2O * 0.5;
	 						}	
						}
						if(proAttenuation == 0)
							exit(1);
						fullIntensity = fullIntensity / proAttenuation;
	 				}
				}
			}
		}
	}
	
	/* If the peptide is 4 residues or INT_4er, look for internal fragment ions. 
	*  (Only +1 right now).  Skip if LCQ data.*/
	fullIntensity = 50;
	if (seqLength > 4 && gParam.fragmentPattern != 'L') {
		for (i = 1; i < seqLength - 2; i++) {
			for (j = i+1; j < seqLength - 1; j++) 
			{
				if(j <= i + 3)
				{
					fragmentMass = gElementMass[HYDROGEN];
					for (k = i; k <= j; k++) fragmentMass += currScorePtr->peptide[k];
			
					if (fragmentMass < msms.scanMassHigh && fragmentMass < parent) 
					{
						if(fragmentMass > lowMassRange && fragmentMass < highMassRange)
						{
							fragmentMass = fragmentMass * 2;
						
							AddPeakToSpectrum(spectrum2, fragmentMass, fullIntensity * INT_FRAG_ATT);
						
							fragmentMass = fragmentMass * 0.5;
						}
					}
				}
			}
		}
	}
	
/*	Add the low mass immonium ions.*/
	for (i = 0; i < seqLength; i++) 
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(gNomMass[j] == currScorePtr->peptide[i])
			{
				if(lowMassIonMass[j] > lowMassRange)
				{
					lowMassIonMass[j] = lowMassIonMass[j] * 2;
					
					AddPeakToSpectrum(spectrum2, lowMassIonMass[j], fullIntensity * lowMassIonIntFactor[j]);
					
					lowMassIonMass[j] = lowMassIonMass[j] * 0.5;
				}
				
				break;
			}
		}
	}
	

/*	Wipe out region where the precursor and derivatives are located; these don't count.*/
	parent    = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	parentMinH2O = (gParam.peptideMW - H2O + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	parentMinNH3 = (gParam.peptideMW - NH3 + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	parentMin2H2O = (gParam.peptideMW - H2O - H2O + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	parentMin2NH3 = (gParam.peptideMW - NH3 - NH3 + (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
	
	parent = parent * 2;
	parentMinH2O = parentMinH2O * 2;
	parentMinNH3 = parentMinNH3 * 2;
	parentMin2H2O = parentMin2H2O * 2;
	parentMin2NH3 = parentMin2NH3 * 2;
	
	spectrum2[((INT_4)(parent + 0.5))] = 0;
	spectrum2[((INT_4)(parent + 0.5)) - 1] = 0;
	spectrum2[((INT_4)(parent + 0.5)) + 1] = 0;
	spectrum2[((INT_4)(parentMinH2O + 0.5))] = 0;
	spectrum2[((INT_4)(parentMinH2O + 0.5)) - 1] = 0;
	spectrum2[((INT_4)(parentMinH2O + 0.5)) + 1] = 0;
	spectrum2[((INT_4)(parentMinNH3 + 0.5))] = 0;
	spectrum2[((INT_4)(parentMinNH3 + 0.5)) - 1] = 0;
	spectrum2[((INT_4)(parentMinNH3 + 0.5)) + 1] = 0;
	spectrum2[((INT_4)(parentMin2H2O + 0.5))] = 0;
	spectrum2[((INT_4)(parentMin2H2O + 0.5)) - 1] = 0;
	spectrum2[((INT_4)(parentMin2H2O + 0.5)) + 1] = 0;
	spectrum2[((INT_4)(parentMin2NH3 + 0.5))] = 0;
	spectrum2[((INT_4)(parentMin2NH3 + 0.5)) - 1] = 0;
	spectrum2[((INT_4)(parentMin2NH3 + 0.5)) + 1] = 0;
	
	if(widePeak)
	{
		spectrum2[((INT_4)(parent + 0.5)) - 2] = 0;
		spectrum2[((INT_4)(parent + 0.5)) + 2] = 0;
		spectrum2[((INT_4)(parentMinH2O + 0.5)) - 2] = 0;
		spectrum2[((INT_4)(parentMinH2O + 0.5)) + 2] = 0;
		spectrum2[((INT_4)(parentMinNH3 + 0.5)) - 2] = 0;
		spectrum2[((INT_4)(parentMinNH3 + 0.5)) + 2] = 0;
		spectrum2[((INT_4)(parentMin2H2O + 0.5)) - 2] = 0;
		spectrum2[((INT_4)(parentMin2H2O + 0.5)) + 2] = 0;
		spectrum2[((INT_4)(parentMin2NH3 + 0.5)) - 2] = 0;
		spectrum2[((INT_4)(parentMin2NH3 + 0.5)) + 2] = 0;
	}
	
/*	Wipe out regions outside of scan range.*/
	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		if(i < ((INT_4)msms.scanMassLow)*2 - 1)
		{
			spectrum2[i] = 0;
		}
		if(i > ((INT_4)msms.scanMassHigh)*2 + 1)
		{
			spectrum2[i] = 0;
		}
	}
	
	/*debug spectrum2*/
	/*if(currScorePtr->rank == 2)
	{
		printf("sequence spectrum \n");	
		for(i = 0; i < SIZEOF_SPECTRA; i++)	
		{
			if(spectrum2[i] != 0)
			{
				printf("%f  %f \n",(REAL_4)i/2,spectrum2[i]);
			}
		}
	}	*/
	
	
	/* Cross-correlation analysis */	
	CrossCorrelate(spectrum2-1, spectrum1-1,(UINT_4) SIZEOF_SPECTRA, tau-1);
		
	/* The cross-correlation score is tau[0] minus the mean of -75 < tau < 75.
	   tau[-1 to -75] are stored in wrapped around order at the end of tau.  */
/*	memcpy(&tau[76], &tau[SIZEOF_SPECTRA - 75], 75 * sizeof(REAL_4));
	intensityAccountedFor = 0.0;
	for (i = 0; i < 150; i++) {
		intensityAccountedFor += tau[i];
	}	
	currScorePtr->crossDressingScore = tau[0] - intensityAccountedFor/150;*/
	
	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		if(tau[i] < 1)
		{
			tau[i] = 0;
		}
	}
	
	/*Since exact matches are exactly symmetrical, I no longer subtract out the average tau from
	-75 to 75, but instead subtract the sum of the differences between points that differ by a factor
	of -1.*/
	intensityAccountedFor = 0;
	for(i = 1; i < 250; i++)	/*compare +1/-1 up to +250/-250*/
	{
		tauDiff = tau[i] - tau[SIZEOF_SPECTRA - i];	/*tau(-250) minus tau(250), etc*/
		if(tauDiff < 0)
		{
			tauDiff = tauDiff * -1;	/*absolute value*/
		}
		intensityAccountedFor += tauDiff;	/*add it all up*/
	}
	currScorePtr->crossDressingScore = tau[0] - intensityAccountedFor/250;/*divide by 250*/
	
	if (spectrum2)	free(spectrum2);	spectrum2 = NULL;	
}
/**************************	AddPeakToSpectrum *******************************************
*
*  
*/
void AddPeakToSpectrum( REAL_4 *spectrum, REAL_4 mass, REAL_4 intensity) {

	char widePeak;
	
	if(intensity < 2)
		return;	/*don't sweat the small stuff*/
		
	if(gParam.fragmentErr <= 0.75)
	{
		widePeak = FALSE;	/*mock peak widths are 1.5 dalton*/
	}
	else
	{
		widePeak = TRUE;	/*mock peak widths are 2.5 daltons*/
	}
	
	/* Make sure that the mass is within the spectrum's range */
	if ( (INT_4)(mass + 0.5) > 2 && (INT_4)(mass + 0.5) < SIZEOF_SPECTRA - 2 ) {
		/* Set the intensity of the peak center */
		if (spectrum[(INT_4) (mass + 0.5)] < intensity)
			spectrum[(INT_4) (mass + 0.5)] = intensity;
			
		/* Set the intensity of the peak sides */
		if (spectrum[((INT_4) (mass + 0.5)) - 1] < intensity * gSidePeakAtt)
			spectrum[((INT_4) (mass + 0.5)) - 1] = intensity * gSidePeakAtt;
		if (spectrum[((INT_4) (mass + 0.5)) + 1] < intensity * gSidePeakAtt)
			spectrum[((INT_4) (mass + 0.5)) + 1] = intensity * gSidePeakAtt;

		if(widePeak)
		{
			/* Set the intensity of the peak side's side*/
			if (spectrum[((INT_4) (mass + 0.5)) - 2] < intensity * gSidePeakAtt * gSidePeakAtt)
				spectrum[((INT_4) (mass + 0.5)) - 2] = intensity * gSidePeakAtt * gSidePeakAtt;
			if (spectrum[((INT_4) (mass + 0.5)) + 2] < intensity * gSidePeakAtt * gSidePeakAtt)
				spectrum[((INT_4) (mass + 0.5)) + 2] = intensity * gSidePeakAtt * gSidePeakAtt;
		}
	}
}		
