/*********************************************************************************************
*  Copyright © 1995-1999      
*  Richard S. Johnson
*  Immunex Corp.
*  Seattle, WA
*   
*  First version: 10/95    
*********************************************************************************************/

/* ANSI Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Lutefisk Headers */
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"

/*	Globals for this file*/
REAL_4	bIonProb				= 0.7;	/*b ion probability		.6	*/
REAL_4	bMinWaterProb			= 0.3;	/*b-18 probability		.3	*/
REAL_4	bMinAmmoniaProb			= 0.15;	/*b-17 probability		.15	*/
REAL_4	bDoublyProbMultiplier	= 0.5;	/*Multiply this for more than one charge		.5	*/
REAL_4	aIonProb				= 0.1;	/*a ion probability		.1	*/
REAL_4	yIonProb				= 0.8;	/*y ion probability		.8	*/
REAL_4	yMinWaterProb			= 0.1;	/*y-18 probability		.1	*/
REAL_4	yMinAmmoniaProb			= 0.1;	/*y-17 probability		.1	*/
REAL_4	yDoublyProbMultiplier	= 0.5;	/*Multiply this for more than one charge		.5	*/
REAL_4	immoniumProb			= 0.2;	/*Immonium ion probability		.2	*/
REAL_4	internalProb			= 0.1;	/*Internal ion probability		.1	*/
REAL_4	internalProProb			= 0.2;	/*Internal ions with N-terminal proline probability		.2	*/


/********************************SequenceScorer*****************************************************
*
*	Assign probability scores to sequences.
*
*/
REAL_4		LutefiskProbScorer(INT_4 *sequence, INT_4 seqLength, REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ)
{
	INT_4 	i, j;
	REAL_4	*randomProb;
	REAL_4	probScore = 0;

/*	Make some space*/
	randomProb = malloc(MAX_ION_NUM * sizeof(REAL_4));
	if(randomProb == NULL)
	{
		printf("SequenceScorer:  Out of memory");
		exit(1);
	}
	/*Initialize*/
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		randomProb[i] 	= 0;
	}
			
/*	Calculate random probability for each ion*/

	CalcRandomProb(randomProb, fragMOverZ, fragNum);
		
/*	Score the sequences*/

	/*initialize ionFound and sequenceProb to all zero, and fill in the precursor related ions*/
//	InitIonFound(ionFound, fragMOverZ, fragNum, randomProb);
	
	/*Get initial probability based on terminal group (Lys and Arg are good; others are not)*/
//	probScore = InitProbScore(sequence, seqLength);
	
//	probScore = FindBIons(ionFound, fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength);
	
	/*Find the y ions*/
//	probScore = FindYIons(ionFound, fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength);
	
	/*Find the internal fragment ions*/
//	probScore = FindInternalIons(ionFound, fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength);
	
	/*Find the immonium ions*/
/*	probScore = FindImmoniumIons(ionFound, fragMOverZ, fragNum, probScore, randomProb,
									sequence, seqLength);*/
		
	/*Change probability score to log base 10 scale*/
	if(probScore > 1)
	{
		probScore = log10(probScore);
	}
	else	/*keep things positive by only logging things over a value of 1*/
	{
		probScore = 0;
	}
	
		
/*	Free array*/
	free(randomProb);
	
	return(probScore);
}


/***********************************FindImmoniumIons******************************************
*
*	Finds and scores the amino acid immonium ions.
*/
REAL_4	FindImmoniumIons(REAL_4 *ionFound, INT_4 *mass, INT_4 ionCount, 
							REAL_4 score, REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength)
{
	REAL_4 lowMassIons[AMINO_ACID_NUMBER][3] = {
	 /* A */     44.0500,     0,     0, 
	 /* R */	 70.0657,  	87.0922, 112.0875,  
	 /* N */	 87.0558,     0,     0,  
	 /* D */	 88.0399,     0,     0, 
	 /* C */		0, 		 0, 	0, 
	 /* E */	102.0555,     0,     0,  
	 /* Q */	 84.0450, 101.0715, 129.0664,
	 /* G */	    0,       0,     0,
	 /* H */	110.0718,     0,     0,  
	 /* I */	 86.0970, 120.0483,     0,  /*I position also represent oxidized Met for qtof*/
	 /* L */	 86.0970,     0,     0,
	 /* K */     84.0814, 101.1079, 129.1028, 
	 /* M */	104.0534,     0,     0, 
	 /* F */	120.0813,     0,     0,  
	 /* P */	 70.0657,     0,     0,  
	 /* S */	 60.0449,     0,     0,
	 /* T */	 74.0606,     0,     0, 
	 /* W */	159.0922,     0,     0, 
	 /* Y */	136.0762,     0,     0,  
	 /* V */	 72.0813,     0,     0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0
	};
	INT_4 i, j, k, amAcidIndex, immoniumIndex;
	REAL_4 individualProb;

	
	for(i = 0; i < seqLength; i++)
	{
		amAcidIndex = -1;
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] <= gMonoMass[j] + gParam.fragmentErr &&
				sequence[i] >= gMonoMass[j] - gParam.fragmentErr)
			{
				amAcidIndex = j;
				break;
			}
		}
		if(amAcidIndex > -1)
		{
			for(j = 0; j < 4; j++)
			{
				if(lowMassIons[amAcidIndex][j] > 0 && lowMassIons[amAcidIndex][j] > mass[0])
				{
					immoniumIndex = 0;
					for(k = 0; k < ionCount; k++)
					{
						if(mass[k] > 160)
						{
							break;
						}
						if(mass[k] < lowMassIons[amAcidIndex][j] + gParam.fragmentErr &&
							mass[k] > lowMassIons[amAcidIndex][j] - gParam.fragmentErr)
						{
							ionFound[k] = 1;
							immoniumIndex = k;
						}
					}
					
					if(immoniumIndex > 0)	/*something was found*/
					{
						individualProb = immoniumProb / randomProb[immoniumIndex];
						score *= individualProb;
					}
					else
					{
						individualProb = (1 - immoniumProb) / (1 - randomProb[immoniumIndex]);
						score *= individualProb;
					}
				}
			}
		}
	}

	return(score);
}
/***********************************FindInternalIons******************************************
*
*	Finds and scores the internal fragment ions.
*/
REAL_4	FindInternalIons(REAL_4 *ionFound, INT_4 *mass, INT_4 ionCount, 
						REAL_4 score, REAL_4 *randomProb, INT_4 *ionRank, INT_4 *sequence, INT_4 seqLength)
{			
	INT_4 i, j, k, residueCount;
	REAL_4 testMass, testFound, individualProb;
	REAL_4 precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN]) / gParam.chargeState;
	BOOLEAN nTermPro, intFragTest;
	
	if(seqLength < 4)
	{
		return(score);	/*need at least four residues for an internal fragment*/
	}
	
	for(i = 1; i < seqLength - 2; i++)
	{
		testMass = sequence[i] + gElementMass[HYDROGEN];
		residueCount = 1;
		if(sequence[i] > gMonoMass[P] - gParam.fragmentErr &&
			sequence[i] < gMonoMass[P] + gParam.fragmentErr)
		{
			nTermPro = TRUE;	/*The N-terminus of this fragment is proline*/
		}
		else
		{
			nTermPro = FALSE;
		}
		for(j = i + 1; j < seqLength - 1; j++)
		{
			testMass += sequence[j];
			residueCount++;
			if(testMass < precursor - gParam.fragmentErr && residueCount < 6
					&& testMass > mass[0])	/*dont bother w/ high mass internal frags*/
			{
				intFragTest = FALSE;
				for(k = 0; k < ionCount; k++)
				{	
					if(mass[k] > testMass + gParam.fragmentErr)
					{
						break;	/*I need to save the k value at the point where this occurs*/
					}
					if(testMass < mass[k] + gParam.fragmentErr &&
						testMass > mass[k] - gParam.fragmentErr)
					{
						if(!nTermPro)
						{
							if((REAL_4)ionRank[k] / (REAL_4)ionCount > 0.5)	/*bottom half of the ions, intensity-wise*/
							{
								intFragTest = TRUE;
								testFound = (REAL_4)ionRank[k] / (REAL_4)ionCount - 0.5; /*ranges from 0 to 0.5*/
								if(testFound > ionFound[k])
								{
									ionFound[k] = testFound;
								}
							}
						}
						else
						{
							intFragTest = TRUE;
							ionFound[k] = 1;	/*internal frags w/ N-terminal Pro are common*/
						}
						break;
					}
				}
				
				/*need to make sure k index is in range*/
				if(k >= ionCount)	
				{
					k = ionCount;
				}
				if(k < 0)
				{
					k = 1;
				}
				
				/*score the probability*/
				if(intFragTest)
				{
					if(nTermPro)
					{
						individualProb = internalProProb / randomProb[k];
						score *= individualProb;
					}
					else
					{
						individualProb = internalProb / randomProb[k];
						score *= individualProb;
					}
				}
				else	/*didn't find any*/
				{
					if(nTermPro)
					{
						individualProb = (1 - internalProProb) / (1 - randomProb[k]);
						score *= individualProb;
					}
					else
					{
						individualProb = (1 - internalProb) / (1 - randomProb[k]);
						score *= individualProb;
					}
				}
			}
		}
	}
	
	return(score);
}

/***********************************InitProbScore**********************************************
*
*	If the C-terminus is Arg or Lys, then give higher probability.
*
*/
REAL_4	InitProbScore(INT_4 *sequence, INT_4 seqLength)
{
	REAL_4 score = 0.05;
	REAL_4 residueMass, testMass;
	INT_4 i;
	BOOLEAN test = FALSE;
	BOOLEAN weirdMass;
	
	
	
	residueMass = sequence[seqLength - 1];
		
	if(residueMass < gMonoMass[R] + gParam.fragmentErr &&
		residueMass > gMonoMass[R] - gParam.fragmentErr)
	{
		score = 0.95;
	}
	else if(residueMass < gMonoMass[K] + gParam.fragmentErr &&
		residueMass > gMonoMass[K] - gParam.fragmentErr)
	{
		score = 0.95;
	}
	else
	{
		for(i = 0; i < gAminoAcidNumber; i++)
		{
			testMass = residueMass - gMonoMass[i];
			if(testMass < gMonoMass[R] + gParam.fragmentErr &&
				testMass > gMonoMass[R] - gParam.fragmentErr)
			{
				score = 0.95;
				break;
			}
			else if(testMass < gMonoMass[K] + gParam.fragmentErr &&
				testMass > gMonoMass[K] - gParam.fragmentErr)
			{
				score = 0.95;
				break;
			}
		}
	}

	return(score);
}

/************************************CalcRandomProb*********************************************
*
*	For each ion, a 200 u window is identified (usually +/- 100 u surrounding it) and the number
*	of ions is counted within the window.  That counted number is divided by the number of possible
*	ions that could fit in that 200 u window, which depends on the instrument resolution.
*/
void	 CalcRandomProb(REAL_4 *randomProb, INT_4 *mass, INT_4 ionCount)
{
	INT_4 i, j, windowCount;
	REAL_4 lowMass, highMass;
	
/*	Initialize*/
	lowMass = mass[0];
	highMass = mass[ionCount-1];
	
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		randomProb[i] = 0;
	}
	
	for(i = 0; i < ionCount; i++)
	{
		windowCount = 0;
		if(mass[i] < lowMass + 100)	/*bottom 200 u window before it moves*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] < lowMass + 200)
				{
					windowCount++;
				}
			}
		}
		else if(mass[i] > highMass - 100)	/*top 200 u window that stops moving*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] > highMass - 200)
				{
					windowCount++;
				}
			}
		}
		else	/*this is the moving window*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] > mass[i] - 100 && mass[j] < mass[i] + 100)
				{
					windowCount++;
				}
			}
		}
		/*calculate the randomness of this ion*/
		randomProb[i] = (REAL_4) windowCount / 200;	/*assuming low resolution of one peak per amu*/
		
	}
	
/*Verify that randomProb nevers equals zero or one (avoid divide by zero later on)*/
	for(i = 0; i < ionCount; i++)
	{
		if(randomProb[i] < 0.005)
		{
			randomProb[i] = 0.005;	/*this is 1 out of 200*/
		}
		if(randomProb[i] > 0.995)
		{
			randomProb[i] = 0.995;	/*this is 199 out of 200*/
		}
	}
	return;
}
/*************************************InitIonFound***********************************************
*
*	Initialize ionFound to zero for each sequence, and set precursor to "found".
*/

void	InitIonFound(REAL_4 *ionFound, INT_4 *mass, INT_4 ionCount, REAL_4 *randomProb)
{
	INT_4 i;
	REAL_4 testMass, water, ammonia;
	
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		ionFound[i]	= 0;
	}
/*	Find precursor ion*/
	testMass = (gParam.peptideMW + gElementMass[HYDROGEN] * gParam.chargeState) / gParam.chargeState;
	water = gElementMass[HYDROGEN] * 2 + gElementMass[OXYGEN];
	ammonia = gElementMass[NITROGEN] + gElementMass[HYDROGEN] * 3;
	for(i = 0; i < ionCount; i++)
	{
		if(mass[i] > testMass + gParam.fragmentErr)
		{
			break;
		}
		if(mass[i] > testMass - water - gParam.fragmentErr)
		{
			if(mass[i] < testMass - water + gParam.fragmentErr)
			{
				ionFound[i] = 1;
				randomProb[i] = 0.005;
			}
			if(mass[i] > testMass - ammonia - gParam.fragmentErr &&
				mass[i] < testMass - ammonia + gParam.fragmentErr)
			{
				ionFound[i] = 1;
				randomProb[i] = 0.005;
			}
			if(mass[i] > testMass - gParam.fragmentErr &&
				mass[i] < testMass + gParam.fragmentErr)
			{
				ionFound[i] = 1;
				randomProb[i] = 0.005;
			}
		}
	}
	return;
}

/**************************************FindBIons*********************************************
*
*	Find the b ions for the sequence and change the ionFound to 1.  Return a value that corresponds
*	to the number of consecutive b ions.
*/

REAL_4 FindBIons(REAL_4 *ionFound, INT_4 *mass, INT_4 ionCount, INT_4 sequenceNum, REAL_4 probScore,
					REAL_4 *randomProb, REAL_4 *intensity, INT_4 *sequence, INT_4 seqLength)
{
	INT_4	ionsInARow = 0;
	INT_4	i, j, k, bIonIndex;
	REAL_4	water, ammonia, bIonTemplate, bIonMin17Template, bIonMin18Template, testFound;
	REAL_4	bIon, bIonMin17, bIonMin18, individualProb, aIon, aIonTemplate, carbonMonoxide;
	BOOLEAN	bIonTest, bMin18Test, bMin17Test, aIonTest, isItAGap;
	
	
/*	Initialize*/
	water = gElementMass[OXYGEN] + gElementMass[HYDROGEN] * 2;
	ammonia = gElementMass[NITROGEN] + gElementMass[HYDROGEN] * 3;
	carbonMonoxide = gElementMass[CARBON] + gElementMass[OXYGEN];
	bIonTemplate = gParam.modifiedNTerm;
	
/*	Start the calculations and searches*/
	for(i = 0; i < seqLength; i++)
	{
		isItAGap = TRUE;
		if(i == 0 || i == seqLength - 1)
		{
			isItAGap = FALSE;	/*won't call the two ends "gaps"*/
		}
		if(isItAGap)
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(sequence[i] < gMonoMass[j] + gParam.fragmentErr &&
					sequence[i] > gMonoMass[j] - gParam.fragmentErr)
				{
					isItAGap = FALSE;
					break;
				}
			}
		}
	
		bIonTemplate += sequence[i];
		bIonMin17Template = bIonTemplate - ammonia;
		bIonMin18Template = bIonTemplate - water;
		aIonTemplate = bIonTemplate - carbonMonoxide;
		for(j = 1; j <= gParam.chargeState; j++)	/*check different charge states*/
		{
			bIon 		= (bIonTemplate + (j-1)*gElementMass[HYDROGEN]) / j;
			bIonMin17 	= (bIonMin17Template + (j-1)*gElementMass[HYDROGEN]) / j;
			bIonMin18 	= (bIonMin18Template + (j-1)*gElementMass[HYDROGEN]) / j;
			aIon		= (aIonTemplate + (j-1)*gElementMass[HYDROGEN]) / j;
			bIonTest = FALSE;
			aIonTest = FALSE;
			bMin18Test = FALSE;
			bMin17Test = FALSE;
			/*apply constraints to charge and mass*/
			if((bIon * j) > ((j-1) * 350) && bIon > mass[0])
			{
				for(k = 0; k < ionCount; k++)
				{
					if(mass[k] > bIon + gParam.fragmentErr)
					{
						break;	/*don't waste any more time looking*/
					}
					if(mass[k] > bIon - gParam.fragmentErr)
					{
						ionFound[k] = 1;
						bIonTest = TRUE;
						bIonIndex = k;
					}
				}
				if(bIonTest)	/*there's a b ion, so look for the b-17 and b-18*/
				{
					k = bIonIndex;
					k--;
					while(mass[k] > aIon - gParam.fragmentErr && k >= 0)
					{
						if(mass[k] > bIonMin17 - gParam.fragmentErr &&
							mass[k] < bIonMin17 + gParam.fragmentErr)
						{
							testFound = intensity[bIonIndex] / (intensity[bIonIndex] + intensity[k]);
							bMin17Test = TRUE;
							if(testFound > ionFound[k])
							{
								ionFound[k] = testFound;
							}
						}
						if(mass[k] > bIonMin18 - gParam.fragmentErr &&
							mass[k] < bIonMin18 + gParam.fragmentErr)
						{
							testFound = intensity[bIonIndex] / (intensity[bIonIndex] + intensity[k]);
							bMin18Test = TRUE;
							if(testFound > ionFound[k])
							{
								ionFound[k] = testFound;
							}
						}
						if(mass[k] > aIon - gParam.fragmentErr &&
							mass[k] < aIon + gParam.fragmentErr)
						{
							testFound = intensity[bIonIndex] / (intensity[bIonIndex] + intensity[k]);
							aIonTest = TRUE;
							if(testFound > ionFound[k])
							{
								ionFound[k] = testFound;
							}
						}
						k--;
					}
				}
				
				/*Calculate the probability scores*/
				if(bIonTest)	/*if the calculated b ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = bIonProb / randomProb[bIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = bIonProb * bDoublyProbMultiplier / randomProb[bIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated b ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1-bIonProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - bIonProb * bDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
				
				if(bMin18Test)	/*if the calculated b-18 ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = bMinWaterProb / randomProb[bIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = bMinWaterProb * bDoublyProbMultiplier / randomProb[bIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1 - bMinWaterProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - bMinWaterProb * bDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
					
				if(bMin17Test)	/*if the calculated b-17 ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = bMinAmmoniaProb / randomProb[bIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = bMinAmmoniaProb * bDoublyProbMultiplier / randomProb[bIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated b ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1 - bMinAmmoniaProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - bMinAmmoniaProb * bDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
				
				if(aIonTest)	/*if the calculated a ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = aIonProb / randomProb[bIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = aIonProb * bDoublyProbMultiplier / randomProb[bIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated b ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1 - aIonProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - aIonProb * bDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
				
				if(isItAGap && j == 1)	/*a gap arises from lack of a y ion, so penalize, but only do it once for j=1*/
				{
					individualProb = (1-bIonProb * 0.5) / (1 - randomProb[k - 1]);	/*0.5 reflects fact that many
																					gaps (except for the ends) are 
																					due to proline*/
					probScore *= individualProb;
				}
			}
		}
	}
	return(probScore);	/*ionsInARow is not used for now*/
}

/**************************************FindYIons*********************************************
*
*	Find the y ions for the sequence and change the ionFound to 1.  Return a value that corresponds
*	to the number of consecutive y ions.
*/

REAL_4 FindYIons(REAL_4 *ionFound, INT_4 *mass, INT_4 ionCount, INT_4 sequenceNum, REAL_4 probScore,
					REAL_4 *randomProb, REAL_4 *intensity, INT_4 *sequence, INT_4 seqLength)
{
	INT_4	ionsInARow = 0;
	INT_4	i, j, k, yIonIndex;
	REAL_4	water, ammonia, yIonTemplate, yIonMin17Template, yIonMin18Template;
	REAL_4	yIon, yIonMin17, yIonMin18, individualProb, testFound;
	BOOLEAN	yIonTest, yMin17Test, yMin18Test, isItAGap;
	
	
/*	Initialize*/
	water = gElementMass[OXYGEN] + gElementMass[HYDROGEN] * 2;
	ammonia = gElementMass[NITROGEN] + gElementMass[HYDROGEN] * 3;
	yIonTemplate = gParam.modifiedCTerm + 2 * gElementMass[HYDROGEN];
	
	for(i = seqLength - 1; i > -1; i--)
	{
		isItAGap = TRUE;
		if(i == 0 || i == seqLength - 1)
		{
			isItAGap = FALSE;	/*won't call the two ends "gaps"*/
		}
		if(isItAGap)
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(sequence[i] < gMonoMass[j] + gParam.fragmentErr &&
					sequence[i] > gMonoMass[j] - gParam.fragmentErr)
				{
					isItAGap = FALSE;
					break;
				}
			}
		}
		
		yIonTemplate += sequence[i];
		yIonMin17Template = yIonTemplate - ammonia;
		yIonMin18Template = yIonTemplate - water;
		for(j = 1; j <= gParam.chargeState; j++)	/*check different charge states*/
		{
			
			yIon 		= (yIonTemplate + (j-1)*gElementMass[HYDROGEN]) / j;
			yIonMin17 	= (yIonMin17Template + (j-1)*gElementMass[HYDROGEN]) / j;
			yIonMin18 	= (yIonMin18Template + (j-1)*gElementMass[HYDROGEN]) / j;
			yIonTest = FALSE;
			yMin17Test = FALSE;
			yMin18Test = FALSE;
			
			/*apply constraints to charge and mass*/
			if((yIon * j) > ((j-1) * 350) && yIon > mass[0])
			{
				for(k = 0; k < ionCount; k++)
				{
					if(mass[k] > yIon + gParam.fragmentErr)
					{
						break;	/*don't waste any more time looking*/
					}
					if(mass[k] > yIon - gParam.fragmentErr)
					{
						ionFound[k] = 1;
						yIonTest = TRUE;
						yIonIndex = k;
					}
				}
				if(yIonTest)	/*there's a y ion, so look for the y-17 and y-18*/
				{
					k = yIonIndex;
					k--;
					while(mass[k] > yIonMin18 - gParam.fragmentErr && k >= 0)
					{
						if(mass[k] > yIonMin17 - gParam.fragmentErr &&
							mass[k] < yIonMin17 + gParam.fragmentErr)
						{	/*y-17 intensity should be less than the y ion*/
							testFound = intensity[yIonIndex] / (intensity[yIonIndex] + intensity[k]);
							yMin17Test = TRUE;
							if(testFound > ionFound[k])
							{
								ionFound[k] = testFound;
							}
						}
						if(mass[k] > yIonMin18 - gParam.fragmentErr &&
							mass[k] < yIonMin18 + gParam.fragmentErr)
						{
							testFound = intensity[yIonIndex] / (intensity[yIonIndex] + intensity[k]);
							yMin18Test = TRUE;
							if(testFound > ionFound[k])
							{
								ionFound[k] = testFound;
							}
						}
						k--;
					}
				}
				
				
				/*Calculate the probability scores*/
				if(yIonTest)	/*if the calculated y ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = yIonProb / randomProb[yIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = yIonProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated y ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1-yIonProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - yIonProb * yDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}

				if(yMin18Test)	/*if the calculated y-18 ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = yMinWaterProb / randomProb[yIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = yMinWaterProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated y ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1 - yMinWaterProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - yMinWaterProb * yDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
				
				if(yMin17Test)	/*if the calculated y-18 ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = yMinAmmoniaProb / randomProb[yIonIndex];
						probScore *= individualProb;
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = yMinAmmoniaProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						probScore *= individualProb;
					}
				}
				else	/*if the calculated y ion is not present*/
				{
					if(k >= ionCount)	/*need to make sure k index is in range*/
					{
						k = ionCount;
					}
					if(k < 0)
					{
						k = 1;
					}
					if(j == 1)
					{
						individualProb = (1 - yMinAmmoniaProb) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
					else
					{
						individualProb = (1 - yMinAmmoniaProb * yDoublyProbMultiplier) / (1 - randomProb[k - 1]);
						probScore *= individualProb;
					}
				}
				if(isItAGap && j == 1)	/*a gap arises from lack of a y ion, so penalize, but only do it once for j=1*/
				{
					individualProb = (1-yIonProb * 0.5) / (1 - randomProb[k - 1]);	/*0.5 reflects fact that many
																					gaps (except for the ends) are 
																					due to proline*/
					probScore *= individualProb;
				}
			}
		}
	}




	return(probScore);
}





