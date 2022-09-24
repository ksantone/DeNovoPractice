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

/*	Richard S. Johnson
	
	Made 12/95
This is for assigning scores to sequences derived from de novo interepretation of CID data. 
  
Input parameters are pointers to the first Sequence struct (containing the completed sequences 
to be scored), and the first MSData struct (containing the cid data - already offsetted), 
"peptideMW" (the peptide molecular weight), "fragmentErr" (the fragment ion m/z tolerance), 
and "chargeState" (the charge state of the precursor ion).  

This function and its related functions all use amino acid residue, element masses, etc that are
100 x the actual value.  This allows for all values to be INT_4s, and this presumably
speeds up the math and the comparisons (although I don't know this for sure).  This bit
of code was originally used by MADMAE, and has been extensively modified for Lutefisk.

All of the functions within this file are either ScoreSequences or are called by 
ScoreSequences.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"



/*
*	Here's some globals that are specific to this file.  They are two amino acid monoisotopic masses
*	times 10000 for Cys, Arg, His, and Lys.  These get modified at the start of the ScoreSequences
*	function in order to accomodate different alkyl groups on cysteine.
*/
INT_4 gCysPlus[AMINO_ACID_NUMBER] = {
		1740463, 2591103, 2170521, 2180361, 2060184, 2320518, 2310678, 1600307, 2400681, 2160933,
		2160933, 2311042, 2340497, 2500776, 2000620, 1900412, 2040569, 2890885, 2660725, 2020776,
		0, 0, 0, 0, 0
	};

INT_4 gArgPlus[AMINO_ACID_NUMBER] = {
	2271382, 3122022, 2701441, 2711281, 2591103, 2851437, 2841597, 2131226, 2931600, 2691852,
	2691852, 2841961, 2871416, 3031695, 2531539, 2431332, 2571488, 3421804, 3191645, 2551695,
	0, 0, 0, 0, 0
};

INT_4 gHisPlus[AMINO_ACID_NUMBER] = {
	2080960, 2931600, 2511018, 2520859, 2400681, 2661015, 2651175, 1940804, 2741178, 2501430,
	2501430, 2651539, 2680994, 2841273, 2341117, 2240909, 2381066, 3231382, 3001222, 2361273,
	0, 0, 0, 0, 0
};

INT_4 gLysPlus[AMINO_ACID_NUMBER] = {
	1991321, 2841961, 2421379, 2431219, 2311042, 2571376, 2561536, 1851164, 2651539, 2411790,
	2411790, 2561899, 2591355, 2751634, 2251477, 2151270, 2291427, 3141743, 2911583, 2271634,
	0, 0, 0, 0, 0
};

INT_4 gProPlus[AMINO_ACID_NUMBER] = {
	1680899, 2531539, 2110957, 2120797, 2000620, 2260954, 2251114, 1540742, 2341117, 2101368,
	2101368, 2251477, 2280933, 2441212, 1941055, 1840848, 1981005, 2831321, 2601161, 1961212,
	0, 0, 0, 0, 0
};

INT_4 gGlnPlus[AMINO_ACID_NUMBER] = {
	1990957, 2841597, 2421015, 2430855, 2310678, 2571012, 2561172, 1850801, 2651175, 2411427,
	2411427, 2561536, 2590991, 2751270, 2251114, 2150906, 2291063, 3141379, 2911219, 2271270,
	0, 0, 0, 0, 0
};

INT_4 gGluPlus[AMINO_ACID_NUMBER] = {
	2000797, 2851437, 2430855, 2440696, 2320518, 2580852, 2571012, 1860641, 2661015, 2421267, 
	2421267, 2571376, 2600831, 2761110, 2260954, 2160746, 2300903, 3151219, 2921059, 2281110,
	0, 0, 0, 0, 0
};

REAL_4 gToleranceNarrow, gToleranceWide;	/*This is used in the "fuzzy logic" of matching calculated and
						observed ion m/z values.*/
						
char gTrypticCterm = FALSE;	/*TRUE if tryptic or lys-c proteolysis with an ion at 147 or 175.*/

INT_4	gCleavageSiteStringent;  /*used for determining quality of spectrum*/
						
char gAmIHere = FALSE;	/*TRUE if tryptic or lys-c proteolysis with an ion at 147 or 175.*/

INT_4 gRightSequence[50] = {	/*Used for checking if the correct sequence is remaining.*/
	7104, 11503, 13706, 9705, 14707, 11308, 14707, 16003, 11308, 12809,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0
};

char 	gDatabaseSeq[MAX_DATABASE_SEQ_NUM][MAX_PEPTIDE_LENGTH];
INT_4 	gPeptideLength[MAX_DATABASE_SEQ_NUM];
INT_4 	gSeqNum = 0;
REAL_8 	gProbScoreMax;
INT_4 	gGapListDipeptideIndex;

/*INT_4 gRightSequence[50] = {	
	1131, 971, 1471, 971, 870, 991, 570, 570, 1010, 570, 570,
	1962, 1150, 2122, 1561, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0
};*/

/*	Globals for haggis in this file*/
REAL_4	bIonProb				= 0.7;	/*b ion probability		.6	*/
REAL_4	bMinWaterProb			= 0.2;	/*b-18 probability		.3	*/
REAL_4	bMinAmmoniaProb			= 0.15;	/*b-17 probability		.15	*/
REAL_4	bMin64IonProb			= 0.15;	/*b-46 probability for b fragments containing oxMet*/
REAL_4	bDoublyProbMultiplier	= 0.5;	/*Multiply this for more than one charge		.5	*/
REAL_4	aIonProb				= 0.1;	/*a ion probability		.1	*/
REAL_4	yIonProb				= 0.8;	/*y ion probability		.8	*/
REAL_4	yMinWaterProb			= 0.1;	/*y-18 probability		.1	*/
REAL_4	yMinAmmoniaProb			= 0.1;	/*y-17 probability		.1	*/
REAL_4	yMin64IonProb			= 0.2;	/*y-46 probability for y fragments containing oxMet*/
REAL_4	yDoublyProbMultiplier	= 0.5;	/*Multiply this for more than one charge		.5	*/
REAL_4	immoniumProb			= 0.2;	/*Immonium ion probability		.2	*/
REAL_4	internalProb			= 0.1;	/*Internal ion probability		.1	*/
REAL_4	internalProProb			= 0.2;	/*Internal ions with N-terminal proline probability		.2	*/



/*****************************GetDatabaseSeq**********************************************
*
*
*
*/
char *GetDatabaseSeq(INT_4 *peptide, INT_4 peptideLength)
{
	char *string, *p;
	BOOLEAN correctSeqTest;
	INT_4 i, j, correctSeq;
	INT_4 metPheDiff = (gMonoMass_x100[F] - gMonoMass_x100[I]) * 1.5;
	INT_4 lysGlnDiff = (gMonoMass_x100[K] - gMonoMass_x100[Q]) * 1.5;
	
	string = (char *) malloc(50);
	if(!string)
	{
		return(NULL);
	}
	
	p = string;
	
	correctSeq = -1;	/*this is an index value, so -1 flags as unassigned*/
	
	if(gParam.qtofErr * gMultiplier <= metPheDiff && gParam.qtofErr != 0)
	{
		metPheDiff = TRUE;
	}
	else
	{
		metPheDiff = FALSE;
	}
	if(gParam.qtofErr * gMultiplier <= lysGlnDiff && gParam.qtofErr != 0)
	{
		lysGlnDiff = TRUE;
	}
	else
	{
		lysGlnDiff = FALSE;
	}
	
	/*Convert the peptide to the actual sequence (rather than index values of gGap*/
	for(i = 0; i < peptideLength; i++)
	{
		if(peptide[i] < gAminoAcidNumber)	/*if not a dipeptide*/
		{
			peptide[i] = gSingAA[peptide[i]];
		}
		else
		{
			free(string);
			return(NULL);	/*something wrong with this sequence, so skip over it*/
		}
	}
	peptide[peptideLength] = 0;		/*NULL terminate the string*/
	
	/*Now compare peptide w/ the database sequences*/
	for(i = 0; i < gSeqNum; i++)
	{
		correctSeqTest = TRUE;
		if(gPeptideLength[i] == peptideLength)
		{
			for(j = 0; j < gPeptideLength[i]; j++)
			{
				if(peptide[j] != gDatabaseSeq[i][j])
				{
					correctSeqTest = FALSE;	/*not the same amino acid*/
					if(peptide[j] == 'L' && gDatabaseSeq[i][j] == 'I')
					{
						correctSeqTest = TRUE;	/*thats because I/L mixed up, so thats ok*/
					}
					if(metPheDiff)
					{
						if(peptide[j] == 'B' && gDatabaseSeq[i][j] == 'm')
						{
							correctSeqTest = TRUE;	/*B means m here (both represent ox Met*/
						}
					}
					else
					{
						if(peptide[j] == 'F' && gDatabaseSeq[i][j] == 'm')
						{
							correctSeqTest = TRUE;
						}
						if(peptide[j] == 'F' && gDatabaseSeq[i][j] == 'B')
						{
							correctSeqTest = TRUE;
						}
					}
					if(!lysGlnDiff)
					{
						if((peptide[j] == 'Q' || peptide[j] == 'K') && 
							(gDatabaseSeq[i][j] == 'Q' || gDatabaseSeq[i][j] == 'K'))
						{
							correctSeqTest = TRUE;
						}
					}
				}
				if(correctSeqTest == FALSE)
				{
					break;
				}
			}
		}
		else
		{
			correctSeqTest = FALSE;
		}
		if(correctSeqTest)
		{
			correctSeq = i;	/*found the correct sequence, now break out*/
			break;
		}
	}

	if(correctSeq >= 0)
	{
		for(i = 0; i < gPeptideLength[correctSeq]; i++)
		{
			p+= sprintf(p, "%c", gDatabaseSeq[correctSeq][i]);
		}
		
	}
	else
	{
		string = NULL;
	}
	
	p+= sprintf(p, "\0");	/* NULL terminate the string */

	return(string);
}

char *PeptideString(INT_4 *peptide, INT_4 peptideLength)
{

	INT_4 i, j, k;
	char  test;
	char *string;
	char *p;
	char *dipeptide;
	INT_4 dipeptideCounter;
	INT_4 dipeptideMass, sum;
	
	string = (char *) malloc(200);
	if(!string)
	{
		return NULL;
	}
	dipeptide = (char *) malloc(3);
	if(dipeptide == NULL)
		return NULL;
	
	p = string;
	
	for(i = 0; i < peptideLength; i++)
	{
		test = TRUE;
		if(peptide[i] < gAminoAcidNumber)	/*if single amino acid, just print the character*/
		{
			p+= sprintf(p, "%c", gSingAA[peptide[i]]);
		}
		else	/*if dipeptide....*/
		{
			if(gGapList[peptide[i]] /gMultiplier < 1000)	/*if the dipeptide is less than 1000 daltons*/
			{
				dipeptideCounter = 0;	/*count the number of possibilities for dipeptides of this mass*/
	
				for(j = 0; j < gAminoAcidNumber; j++)
				{
					for(k = j; k < gAminoAcidNumber; k++)
					{
						sum = gMonoMass_x100[j] + gMonoMass_x100[k];
						if(sum <= gGapList[peptide[i]] + gToleranceWide &&
							sum >= gGapList[peptide[i]] - gToleranceWide)
						{
							dipeptideCounter++;
						}
					}
				}
				if(dipeptideCounter == 0)
				{
					dipeptideCounter = 1;
				}
				if(dipeptideCounter > 1 || gParam.qtofErr == 0)	/*if there is more than one choice, 
																or if the qtofErr was not used,
																then print the mass in brackets*/
				{
					dipeptideMass = gGapList[peptide[i]];
					/*if(gNodeCorrection[peptide[i]] > 5)
					{
						dipeptideMass = dipeptideMass + 1;
					}
					else if(gNodeCorrection[peptide[i]] <= -5)
					{
						gNodeCorrection[i] = gNodeCorrection[i] + 10;
						dipeptideMass = dipeptideMass - 1;
					}*/
					
					
					if(gMultiplier == 1)
					{
						p+= sprintf(p, "[%3d]", (dipeptideMass)/gMultiplier);
					}
					else if(gMultiplier == 10)
					{
						p+= sprintf(p, "[%.1f]", ((REAL_4)dipeptideMass)/gMultiplier);
					}
					else if(gMultiplier == 100)
					{
						p+= sprintf(p, "[%.2f]", ((REAL_4)dipeptideMass)/gMultiplier);
					}
					else
					{
						p+= sprintf(p, "[%.3f]", ((REAL_4)dipeptideMass)/gMultiplier);
					}
				}
				else	/*if there is only one choice, then print the two amino acids in brackets*/
				{
					test = TRUE;
					for(j = 0; j < gAminoAcidNumber; j++)
					{
						for(k = j; k < gAminoAcidNumber; k++)
						{
							if(gGapList[peptide[i]] == gMonoMass_x100[j] + gMonoMass_x100[k])
							{
								dipeptide[0] = gSingAA[j];
								dipeptide[1] = gSingAA[k];
								dipeptide[2] = 0;
								test = FALSE;
							}
						}
					}
					if(test)	/*non-standard dipeptide*/
					{
						if(gMultiplier == 1)
						{
							p+= sprintf(p, "[%3d]", (gGapList[peptide[i]])/gMultiplier);
						}
						else if(gMultiplier == 10)
						{
							p+= sprintf(p, "[%.1f]", ((REAL_4)gGapList[peptide[i]])/gMultiplier);
						}
						else if(gMultiplier == 100)
						{
							p+= sprintf(p, "[%.2f]", ((REAL_4)gGapList[peptide[i]])/gMultiplier);
						}
						else
						{
							p+= sprintf(p, "[%.3f]", ((REAL_4)gGapList[peptide[i]])/gMultiplier);
						}
					}
					else	/*standard dipeptide*/
					{
						p += sprintf(p, "[%s]", dipeptide);
					}
				}
			}
			else	/*if its some sort of weird dipeptide of mass greater than 1000 Da, then dump it an integer*/
			{
				p+= sprintf(p, "[%3d]", gGapList[peptide[i]] / gMultiplier);
			}
		}
	}
	p+= sprintf(p, "\0");	/* NULL terminate the string */
	
	return string;
}	


/****************************AddPyroGlu**********************************************
*
*
*
*/
void AddPyroGlu()
{
	
	gSingAA[gAminoAcidNumber] = 'q';
	gMonoMass[gAminoAcidNumber] = 111.032028;
	gAvMass[gAminoAcidNumber] = 111.05;
	gNomMass[gAminoAcidNumber] = 111;
	gMonoMass_x100[gAminoAcidNumber] = gMonoMass[gAminoAcidNumber] * gMultiplier + 0.5;
	gAminoAcidNumber++;
	return;	
}

/***************************AddDatabaseSequences*************************************
*
*	This bit of code is for adding sequence(s) derived from database matches using
*	sequest, bequest, mascot, etc.  The idea is to add database-derived sequences
*	to the de novo - derived sequences to have them battle it out in the final 
*	scoring and ranking.  Presumably, if the database sequences are correct, they will
*	best account for the data, compared to the de novo sequences; if they are wrong, then
*	the de novo sequences will win out.
*
*/

void AddDatabaseSequences(struct Sequence *firstSequencePtr)
{
	struct Sequence *currPtr, *newPtr;
	char *stringBuffer, test;
	char databaseSequences[MAX_DATABASE_SEQ_NUM][MAX_PEPTIDE_LENGTH];
	REAL_4 nNodeValue, peptideStartMass, peptideMass;	
	INT_4 peptideLength[MAX_DATABASE_SEQ_NUM];
	INT_4 i, j, k, seqNum, intNodeValue;
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	REAL_4 peptideNodeMass;
	FILE *fp;
	INT_4 lysGlnDiff = (gMonoMass_x100[K] - gMonoMass_x100[Q]) * 1.5;
	
	if(gParam.qtofErr <= lysGlnDiff && gParam.qtofErr != 0)
	{
		lysGlnDiff = TRUE;
	}
	else
	{
		lysGlnDiff = FALSE;
	}
	peptideStartMass = gParam.modifiedCTerm + gParam.modifiedNTerm;
	
	/*Determine the initial nodeValue, depending on N-terminal modification*/
	
	nNodeValue = gParam.modifiedNTerm / gMultiplier;

	stringBuffer = (char *)malloc(258);
    if (stringBuffer == NULL)
    {
        printf("LutefiskScore: stringBuffer == NULL\n");
        exit(1);
    }
    
    fp = fopen(gParam.databaseSequences,"r");
    if (fp == NULL)
    {
   		printf("Cannot open the database sequence file '%s'\n", gParam.databaseSequences);
    	exit(1);
    }
    
    if(gParam.fMonitor)
    {
       	printf("Processing database sequence file '%s'\n", gParam.databaseSequences);
    }
	
	/*Read the information into a character array of single letter amino acid codes.*/
	test = TRUE;
	seqNum = 0;
	gSeqNum = 0;
	while(!feof(fp) && test)
	{
		for(i = 0; i < 258; i++)
		{
			stringBuffer[i] = 0;
		}
		test = FALSE;
		if(my_fgets(stringBuffer, 256, fp) == NULL)
		{
			continue;
		}
		
		j = 0;
		
		/*the gDatabaseSeq and related globals keep track of the actual characters from the
		database sequences w/o conversion of I to L or B to F.  This info is retained in order
		to have the output have the same sequence as in the database.  N-terminal 'q' is considered
		as pyroglutamic acid, but 'q' anywhere else is converted to 'Q' and considered to be Gln*/
		
		while(stringBuffer[j] >= 65 && stringBuffer[j] <= 121)
		{
			gDatabaseSeq[gSeqNum][j] = stringBuffer[j];
			/*Except for m, convert lower case to upper case for global database seqs*/
			if(gDatabaseSeq[gSeqNum][j] >= 97)
			{
				if(gDatabaseSeq[gSeqNum][j] != 109)
				{
					if(gDatabaseSeq[gSeqNum][j] != 113 || j != 0)	/*except for N-terminal q (pyroGlu)*/
					{
						gDatabaseSeq[gSeqNum][j] = gDatabaseSeq[gSeqNum][j] - 32;
					}
				}
			}
			test = TRUE;
			
			/*Now that the original stringBuffer was saved in the global array, I muck w/ it by converting
			m to F and capitalizing everything else except for q*/
			if(stringBuffer[j] >= 97)
			{
				if(stringBuffer[j] == 109)
				{
					stringBuffer[j] = 70;	/*m (109) usually means oxidized Met, which is almost
												the same mass as F (70)*/
				}
				else if(stringBuffer[j] == 113 && j == 0)	/*q (113) means N-terminal pyroglu*/
				{
					AddPyroGlu();	/*adds pyroGlu to list*/
				}
				else
				{
					stringBuffer[j] = stringBuffer[j] - 32;
				}
			}
			/*If tolerance can't differentiate lys from gln, then set Q's to K*/
			if(stringBuffer[j] == 81 && !lysGlnDiff)
			{
				stringBuffer[j] = 75;
			}
			if(stringBuffer[j] == 66) /*Change B to F*/
			{
				stringBuffer[j] = 70;
			}
			if(stringBuffer[j] == 73)
			{
				stringBuffer[j] = 76;
			}
			databaseSequences[seqNum][j] = stringBuffer[j];
			j+=1;
		}
		if(j > 0)	/*a legitimate sequence found*/
		{
			for(i = j; i < 258; i++)
			{
				databaseSequences[seqNum][i] = 0;	/*fill in the rest w/ zero's*/
				gDatabaseSeq[gSeqNum][i] = 0;
			}
			peptideLength[seqNum] = j;
			gPeptideLength[gSeqNum] = j;
			gSeqNum += 1;
			seqNum += 1;
		}
	}
	
	fclose(fp);
	
/*	Convert the characters in databaseSequences to monoisotopic masses*/
	for(i = 0; i < seqNum; i++)
	{
		for(j = 0; j < MAX_PEPTIDE_LENGTH; j++)
		{
			peptide[j] = 0;
		}
		peptideMass = peptideStartMass;
		peptideNodeMass = nNodeValue;
		for(j = 0; j < peptideLength[i]; j++)
		{
			for(k = 0; k < gAminoAcidNumber; k++)
			{
				if(gSingAA[k] == databaseSequences[i][j])
				{
					peptideNodeMass += gMonoMass[k];
					peptide[j] = gMonoMass_x100[k];
					peptideMass += gMonoMass_x100[k];
					break;
				}
			}
		}
	
		peptideNodeMass = peptideNodeMass * gMultiplier + 0.5;
		intNodeValue = peptideNodeMass;
    	if(peptideMass <= gParam.peptideMW + gParam.peptideErr &&
    		peptideMass >= gParam.peptideMW - gParam.peptideErr)
    	{
			/*first find the end of the list of the de novo sequences*/
			currPtr = firstSequencePtr;
			while(currPtr->next != NULL)
			{
				currPtr = currPtr->next;
			}
			
			newPtr = (struct Sequence *) malloc(sizeof(struct Sequence));
			if(newPtr == NULL)
			{
				printf("newPtr == NULL");
				exit(1);
			}
		
			currPtr->next = newPtr;
			newPtr->next = NULL;
			for(j = 0; j < peptideLength[i]; j++)
			{
				newPtr->peptide[j] = peptide[j];
			}
			newPtr->peptideLength = peptideLength[i];
			newPtr->score = 500;
			newPtr->nodeValue = intNodeValue;
			newPtr->nodeCorrection = 0;
			newPtr->gapNum = -100;	/*FLAG FOR BEING A DATABASE SEQUENCE*/
		}
	}
	
	printf("%ld database sequences added.\n", seqNum);
	free(stringBuffer);
	return;
}
/***************************ScoreAttenuationFromCalfactor****************************
*
*	Sequences that need calFactors that deviate from 1.0000000000 are attenuated
*	according to the degree to which they differ.
*/
REAL_4	ScoreAttenuationFromCalfactor(REAL_4 calFactor, REAL_4 intScore)
{
	REAL_4 newIntScore = 1;
	REAL_4 fudgeFactor;
	
	if(calFactor <= 0)	/*calFactor = 0 was a flag for the output that no recalibration
						was done to the data due to a lack of sufficient calibrant ions.*/
		return(intScore);
	if(calFactor > 1)
	{
		fudgeFactor = 1 / calFactor;
	}
	else
	{
		fudgeFactor = calFactor;
	}
	fudgeFactor = (fudgeFactor - 0.99) * 100;
	newIntScore = intScore * fudgeFactor;
	

	return(newIntScore);
}
/******************************ProlineInternalFrag************************************************
*
*	ProlineInternalFrag identifies internal fragment ions that contain proline at the
*	N-terminus.  
*/
void ProlineInternalFrag(REAL_4 *ionFound, INT_4 *fragMOverZ, 
				INT_4 *sequence, INT_4 seqLength, INT_4 fragNum)
{
	INT_4 i, j, k, intFrag, intFragMinErr, intFragPlusErr;
	INT_4 massDiff, precursor;
	REAL_4 currentIonFound;
	
	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) 
				/ gParam.chargeState;

	if(seqLength >= 4)	/*Sequences less than 4 amino acids cannot have internal fragment ions.*/
	{
		for(i = 1; i < (seqLength - 3); i++)	/*This is the N-terminus of the fragment.*/
		{
			if(sequence[i] == gGapList[P])
			{
				for(j = (i + 1); j < (seqLength - 2); j++)	/*C-terminus of the fragment.*/
				{
					if(j < i + 3)	/*just look at int frags less than 4 aa*/
					{
						intFrag = gElementMass_x100[HYDROGEN];	/*Calc the mass of the 
																int frag.*/
						for(k = i; k <= j; k++)
						{
							intFrag += sequence[k];
						}
						intFragMinErr = intFrag - gToleranceWide;
						intFragPlusErr = intFrag + gToleranceWide;
				
						if(intFragMinErr < precursor)	/*Only count those matches where 
														the internal frag mass is less 
														than the precursor m/z value.*/
						{
							k = 0;
							/*Look for this ion.*/
							while(fragMOverZ[k] <= intFragPlusErr && k < fragNum)
							{
								if(fragMOverZ[k] >= intFragMinErr)
								{
									currentIonFound = ionFound[k];
									massDiff = abs(intFrag - fragMOverZ[k]);
									ionFound[k] = CalcIonFound(ionFound[k], massDiff);
									/*Attenuate the internal ion intensity.*/
									ionFound[k] = ionFound[k] * INTERNAL_FRAG_MULTIPLIER ;
									if(currentIonFound > ionFound[k])
									{
										ionFound[k] = currentIonFound;
									}
								}
								k++;
							}
						}
					}	/*if(j < i + 3)*/
				}	/*for(j = (i + 1); j < (seqLength - 2); j++)*/
			}	/*if(sequence[i] == gGapList[P])*/
		}	/*for(i = 1; i < (seqLength - 3); i++)*/
	}	/*if(seqLength >= 4)*/

	return;	
}

/**************************Recalibrate************************************************
*
*	If qtof data is used and if the qtofErr feature employed, then the mass data 
*	is recalibrated according to the errors found between the calculated y and b
*	ions and the original data.  The original mass values are restored at the end
*	of the scoring loop.
*
*/

REAL_4	Recalibrate(INT_4 fragNum, INT_4 *fragMOverZ, INT_4 *sequence, INT_4 seqLength,
					INT_4 *fragIntensity)
{
	INT_4 i, j, k, bCal, yCal, errNum;
	INT_4 bCalStart, yCalStart;
	INT_4 bIonMass, bIonMassMinErr, bIonMassPlusErr;
	INT_4 yIonMass, yIonMassMinErr, yIonMassPlusErr;
	INT_4 *ionFound;
	INT_4 yCalCorrection = 0;
	INT_4 bCalCorrection = 0;
	INT_4 chargeLimit = 1;
	REAL_8 *byError, stDev, avCorrectionFactor;
	REAL_4 precursor, avIntensity;
	REAL_4 totalIntensity = 0;
	REAL_4 totalErrorIntensity = 0;
	REAL_4 ratio, offSet;
	char test;
			
/*	Initialize.*/
	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) /
				gParam.chargeState;
				
	/*For +3 precursors, only look at +1 and +2 fragment ions.*/
	if(gParam.chargeState > 1)
	{
		chargeLimit = gParam.chargeState - 1;
	}
	
	stDev = 0;
	
	ionFound = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(ionFound == NULL)
	{
		printf("Recalibrate: out of memory.\n");
		exit(1);
	}
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		ionFound[i] = 0;
	}
	byError = (double *) malloc(MAX_ION_NUM * sizeof(REAL_8));
	if(byError == NULL)
	{
		printf("Recalibrate: out of memory.\n");
		exit(1);
	}
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		byError[i] = 0;
	}
	avIntensity = 0;
	for(i = 0; i < fragNum; i++)
	{
		avIntensity += fragIntensity[i];
		totalIntensity += fragIntensity[i];
	}
	if(fragNum == 0)
	{
		printf("Recalibrate: fragNum = 0\n");
		exit(1);
	}
	avIntensity = avIntensity / fragNum;
	

/*	
	Initialize the starting b ion mass (acetylated, etc).
*/
	bCalStart = gParam.modifiedNTerm + 0.5;
	/*Determine the correction factor for high accuracy*/
	i = gParam.modifiedNTerm * 10 + 0.5;
	bCalCorrection = i - bCalStart * 10;

/*	
	Initialize the starting mass for y ions.
*/
	yCalStart = (2 * gElementMass_x100[HYDROGEN]) + gParam.modifiedCTerm + 0.5;
	i = (gParam.modifiedCTerm + (2 * gElementMass[HYDROGEN] * gMultiplier)) * 10 + 0.5;
	yCalCorrection = i - yCalStart * 10;
	

	for(i = (seqLength - 1); i > 0; i--)	/*Don't do this loop for i = 0 (doesnt make sense).*/
	{
/*
	Calculate the singly charged y ion mass.
*/
		yCal = YCalculator(i, sequence, seqLength, yCalStart, yCalCorrection);	
			
/*
	Calculate the singly charged b ion mass.
*/
		bCal = BCalculator(i, sequence, bCalStart, bCalCorrection);		
		
		for(j = 1; j <= chargeLimit; j++)
		{
				
		/*bIonMass = bCal;*/
			bIonMass = (bCal + (j * gElementMass_x100[HYDROGEN]) - 
						gElementMass_x100[HYDROGEN]) / j;
			bIonMassMinErr = bIonMass - gParam.fragmentErr;
			bIonMassPlusErr = bIonMass + gParam.fragmentErr;
		/*yIonMass = yCal;*/
			yIonMass = (yCal + (j * gElementMass_x100[HYDROGEN]) 
						- gElementMass_x100[HYDROGEN]) / j;
			yIonMassMinErr = yIonMass - gParam.fragmentErr;
			yIonMassPlusErr = yIonMass + gParam.fragmentErr;
			
/*	Search for b ions.*/
			if(j == 1) /*Only look for +1 b ions*/
			{
				k = fragNum - 1;
				while(fragMOverZ[k] >= bIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= bIonMassPlusErr)
					{
						if(ionFound[k] == 0)
						{
							ionFound[k] = bIonMass;
						}
						else
						{
							ionFound[k] = 0;	/*if both b and y ions match to same ion*/
						}
					}
					k--;
				}
			}
				
/*	Search for the y ion values.*/
			test = TRUE;	/*Only look for multiply charged y ions greater than precursor*/
			if(j > 1)
			{
				if(yIonMass <= precursor + gParam.fragmentErr)
				{
					test = FALSE;
				}
			}
			if(test)
			{
				k = fragNum - 1;
				while(fragMOverZ[k] >= yIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= yIonMassPlusErr)
					{
						if(ionFound[k] == 0)
						{
							ionFound[k] = yIonMass ;
						}
						else
						{
							ionFound[k] = 0;	/*if both b and y ions match to same ion*/
						}
					}
					k--;
				}
			}	/*if(test)*/
		}	/*for j*/
	}	/*for i*/
	
/*	Determine the calibration slope change*/

/*	First find ions greater than the precursor mass or if below the precursor has unusually
	high intensity.  Don't use ions of unusually low abundance.
*/
	
	avCorrectionFactor = 0;
	errNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(ionFound[i] != 0 && fragMOverZ[i] != 0 && 
			(fragMOverZ[i] > precursor || fragIntensity[i] > 2 * avIntensity)
			&&	fragIntensity[i] > avIntensity / 4)
		{
			ratio = (REAL_4)ionFound[i] / (REAL_4)fragMOverZ[i];
			if(ratio > 0.9995 && ratio < 1.0005)
			{
				byError[i] = ratio;
				errNum++;
				avCorrectionFactor += byError[i];
			}
			else
			{
				byError[i] = 100;
			}
		}
		else
		{
			byError[i] = 100;
		}
	}
	if(errNum < 3)	/*if less than three points for determining correction, then dont
					recalibrate*/
	{
		free(byError);
		free(ionFound);
		return(0);	/*0 is returned, but is not applied to data*/
	}
	
/*	Eliminate any outliers.*/
	stDev = StandardDeviationOfTheBYErrors(byError, fragNum);
	/*stDev = stDev * 2;*/
	avCorrectionFactor = avCorrectionFactor / errNum;
	if(stDev != 0)
	{
		for(i = 0; i < fragNum; i++)
		{
			if(byError[i] != 100)
			{
				if(byError[i] < avCorrectionFactor - stDev || 
						byError[i] > avCorrectionFactor + stDev)
				{
					byError[i] = 100;
				}
			}
		}
	}
	
/*	Determine the correction factor.*/
	avCorrectionFactor = 0;
	errNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(byError[i] < 100)
		{
			avCorrectionFactor += (byError[i] * fragIntensity[i]);
			totalErrorIntensity += fragIntensity[i];
			errNum++;
		}
	}
	if(errNum < 2)	/*if less than two points for determining correction, then dont
					recalibrate*/
	{
		free(byError);
		free(ionFound);
		return(0);	/*0 is returned, but is not applied to data*/
	}
	if(totalErrorIntensity == 0)
	{
		printf("Recalibrate: totalErrorIntensity = 0\n");
		exit(1);
	}
	avCorrectionFactor = avCorrectionFactor / totalErrorIntensity;
	
/*	Don't use correction factors greater than 500 ppm.*/
	if(avCorrectionFactor > 1.0005 || avCorrectionFactor < 0.9995)
	{
		avCorrectionFactor = 1;
	}
	
/*	Apply the correction factor to the masses.*/
	for(i = 0; i < fragNum; i++)
	{
		fragMOverZ[i] = fragMOverZ[i] * avCorrectionFactor;
		if(ionFound[i] != 0)
		{
			byError[i] = ionFound[i] - fragMOverZ[i];
		}
		else
		{
			byError[i] = 100;
		}
	}
	
/*	Determine the slope offset and apply*/
	offSet = 0;
	errNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(byError[i] != 100)
		{
			offSet += byError[i];
			errNum++;
		}
	}
	
	if(offSet > 0)
	{
		offSet = offSet / errNum + 0.5;
	}
	else
	{
		offSet = offSet / errNum - 0.5;
	}
	/*this offset is only for very minor adjustments and limited to +/- 1  unit of mass*/
	if(offSet >= 1)
	{
		offSet = 1;
	}
	else if(offSet <= -1)
	{
		offSet = -1;
	}
	else
	{
		offSet = 0;
	}
	
	for(i = 0; i < fragNum; i++)
	{
		fragMOverZ[i] = fragMOverZ[i] + offSet;
	}
	
	free(ionFound);
	free(byError);
	return(avCorrectionFactor);
}

/**********************************IsThisADuplicate************************************
*
*	Find out if this sequence is similar to existing sequences that have been scored
*	and stored (ie, same intOnlyScore and where elements of the sequence are identical
*	or within tolerance).
*/
char IsThisADuplicate(struct SequenceScore *firstScorePtr, INT_4 *sequence, 
						REAL_4 intOnlyScore, REAL_4 intScore, INT_4 seqLength)
{
	struct SequenceScore *currPtr;
	INT_4 i, j, k, m, n;
	char test = TRUE;
	
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		test = TRUE;
		j = intOnlyScore * 1000 + 0.5;
		k = (currPtr->intensityOnlyScore) * 1000 + 0.5;
		m = intScore * 1000 + 0.5;
		n = (currPtr->intensityScore) * 1000 + 0.5;
		if(j == k && m == n)
		{
			i = 0;
			while(currPtr->peptide[i] != 0)
			{
				i++;
			}
			if(i == seqLength)
			{
				for(i = 0; i < seqLength; i++)
				{
					if(sequence[i] < currPtr->peptide[i] - gToleranceNarrow ||
						sequence[i] > currPtr->peptide[i] + gToleranceNarrow)
					{
						test = FALSE;
					}
					
				}
				if(test)
				{
					test = FALSE;
					return(test);
				}
				
			}
		
		}
	
		currPtr = currPtr->next;
	}


	test = TRUE;	/*if it made it this far, then its not a duplicate*/
	return(test);
}
/***********************************GoodSequence***************************************
*
*	This function determines if the "new" sequence is just a repeat of an old sequence,
*	and it also makes sure that the mass of the sequence equals the peptide mass.
*/

char	GoodSequence(struct Sequence *firstSequencePtr, INT_4 *peptide, INT_4 *aaCorrection, 
					INT_4 peptideLength)
{
	INT_4 peptideMass, peptideMassCorrection, i, j;
	char test = TRUE;	/*test is assumed to be true, unless proven otherwise, and this is the returned value*/
	struct Sequence *checkPtr;
	

	checkPtr = firstSequencePtr;
	while(checkPtr != NULL && test)
	{
		test = FALSE;
		for(i = 0; i < checkPtr->peptideLength; i++)
		{
			if(peptide[i] != checkPtr->peptide[i])
			{
				test = TRUE;
			}
		}
		checkPtr = checkPtr->next;
	}
	if(test)	/*verify that the sequences matches the peptide molecular weight.*/
	{
		
		peptideMass = gParam.modifiedNTerm + 0.5;
		/*Determine the correction factor for high accuracy*/
		i = gParam.modifiedNTerm * 10 + 0.5;
		peptideMassCorrection = i - peptideMass * 10;

		peptideMass += (gParam.modifiedCTerm + 0.5);
		i = gParam.modifiedCTerm * 10 + 0.5;
		j = gParam.modifiedCTerm + 0.5;
		j *= 10;
		peptideMassCorrection += (i - j);
		
		for(i = 0; i < peptideLength; i++)
		{
			peptideMass += peptide[i];
			peptideMassCorrection += aaCorrection[i];
			if(peptideMassCorrection >= 10)
			{
				peptideMass += 1;
				peptideMassCorrection -= 10;
			}
			if(peptideMassCorrection <= -10)
			{
				peptideMass -= 1;
				peptideMassCorrection += 10;
			}
		}
		if(peptideMass < gParam.peptideMW - gParam.peptideErr || 
			peptideMass > gParam.peptideMW + gParam.peptideErr)
		{
			test = FALSE;
		}
	}
	return(test);
}


/*****************************ExpandSequences*******************************************
*
*	This function takes the final list of completed sequences and (for qtof only) expands
*	the number of sequences depending on the other members of the updated gGapList that are
*	within the fragment ion tolerance that is set wide during the subsequence buildup.  Later 
*	these sequences (the original ones plus the new ones based upon the old sequences) are scored
*	using tighter qtof error tolerances and the data is recalibrated for each sequence.
*/

void	ExpandSequences(struct Sequence *firstSequencePtr)
{
	struct Sequence *currPtr = NULL;
	struct Sequence *stopPtr = NULL;
	struct Sequence *newPtr = NULL;
	struct Sequence *lastPtr = NULL;
	struct Sequence *checkPtr = NULL;
	struct Sequence *testPtr = NULL; /*debug*/
	INT_4 i, j, peptideNum;
	INT_4 n = 1, p = 0, q = 0; /*debug*/
	INT_4 sequenceMatrix[MAX_PEPTIDE_LENGTH][30], aaNum[MAX_PEPTIDE_LENGTH];
	INT_4 sequence[MAX_PEPTIDE_LENGTH];
	char stop = TRUE;
	char allNegative;
	char test;
	INT_4 cycle, countTheSeqs;
	INT_4 peptide[MAX_PEPTIDE_LENGTH], aaCorrection[MAX_PEPTIDE_LENGTH];
	INT_4 peptideLength, score, nodeValue, gapNum;
	INT_2 nodeCorrection;
	
	currPtr = firstSequencePtr;
	if(currPtr == NULL)
	{
		printf("There are no sequences in ExpandSequences.\n");
		exit(1);
	}

	while(currPtr->next != NULL)
	{
		currPtr = currPtr->next;
	}
	stopPtr = currPtr;
	lastPtr = currPtr;
	currPtr = firstSequencePtr;
	
	while(stop)
	{
		if(currPtr == stopPtr)	/*Signal that the last of the original sequences is being examined.*/
			stop = FALSE;
			
		/*initialize variables.*/
		for(i = 0; i < MAX_PEPTIDE_LENGTH; i++)
		{
			aaNum[i] = 0;
			sequence[i] = 0;
		}
		for(i = 0; i < MAX_PEPTIDE_LENGTH; i++)
		{
			for(j = 0; j < 30; j++)
			{
				sequenceMatrix[i][j] = 0;
			}
		}
		
/*	Determine the different residue possibilities for each position in the sequence.*/

		for(i = 0; i < currPtr->peptideLength; i++)
		{
			for(j = 0; j <= gGapListIndex; j++)
			{
				if(currPtr->peptide[i] <= gGapList[j] + gParam.fragmentErr &&
				currPtr->peptide[i] >= gGapList[j] - gParam.fragmentErr)
				{
					test = FALSE;
					sequenceMatrix[i][aaNum[i]] = j;
					aaNum[i] = aaNum[i] + 1;
					if(aaNum[i] >= 30)
					{
						printf("LutefiskScore: Number of residues in sequenceMatrix exceeds 30.\n");
						exit(1);
					}
				}
			}
			if(aaNum[i] == 0)
			{
				gGapListIndex++;
				gGapList[gGapListIndex] = currPtr->peptide[i];
				gNodeCorrection[gGapListIndex] = 0;
				sequenceMatrix[i][aaNum[i]] = gGapListIndex;
				aaNum[i] += 1;
			}
			n *= aaNum[i];
		}
		if(gParam.proteolysis == 'T')	/*If the C-terminal amino acid is close to R or K, then force it for tryptics*/
		{
			test = FALSE;
			for(i = 0; i < aaNum[(currPtr->peptideLength) - 1]; i++)
			{
				if(gGapList[sequenceMatrix[(currPtr->peptideLength)-1][i]] <= gGapList[K] + gParam.fragmentErr &&
					gGapList[sequenceMatrix[(currPtr->peptideLength)-1][i]] >= gGapList[K] - gParam.fragmentErr)
				{
					test = TRUE;
				}
			}
			if(test)
			{
				aaNum[(currPtr->peptideLength)-1] = 1;
				sequenceMatrix[(currPtr->peptideLength)-1][0] = K;
			}
			test = FALSE;
			for(i = 0; i < aaNum[(currPtr->peptideLength)-1]; i++)
			{
				if(gGapList[sequenceMatrix[(currPtr->peptideLength)-1][i]] <= gGapList[R] + gParam.fragmentErr &&
					gGapList[sequenceMatrix[(currPtr->peptideLength)-1][i]] >= gGapList[R] - gParam.fragmentErr)
				{
					test = TRUE;
				}
			}
			if(test)
			{
				aaNum[(currPtr->peptideLength)-1] = 1;
				sequenceMatrix[(currPtr->peptideLength)-1][0] = R;
			}
		}

/*	Test that there won't be too many sequences made.  If too many, then eliminate all multiple
	choices for each position, except for Q/K and m/F.*/
		
		peptideNum = 1;
		for(i = 0; i < currPtr->peptideLength; i++)
		{
			peptideNum *= aaNum[i];
		}
		
		if((peptideNum > 250 && gCorrectMass) || (peptideNum > 50 && !gCorrectMass))	/*100 is arbitrary  gAminoAcidNumber*/
		{
			/*anything thats not a single aa is made negative*/
			for(i = 0; i < currPtr->peptideLength; i++)
			{
				if(aaNum[i] > 1)
				{
					for(j = 0; j < aaNum[i]; j++)
					{
						if(sequenceMatrix[i][j] >= gAminoAcidNumber)
						{
							sequenceMatrix[i][j] = -1 * sequenceMatrix[i][j];
						}
					}
				}
			}
			/*if all are negative then the first is made positive*/
			for(i = 0; i < currPtr->peptideLength; i++)
			{
				if(aaNum[i] > 1)
				{
					allNegative = TRUE;
					for(j = 0; j < aaNum[i]; j++)
					{
						if(sequenceMatrix[i][j] > 0)
						{
							allNegative = FALSE;
						}
					}
					if(allNegative)
					{
						sequenceMatrix[i][0] = -1 * sequenceMatrix[i][0];
					}
				}
			}
			/*eliminate negatives*/
 			for(i = 0; i < currPtr->peptideLength; i++)
			{
				for(j = 0; j < aaNum[i]; j++)
				{
					if(sequenceMatrix[i][j] < 0)
						break;
				}
				aaNum[i] = j;
			}
							
		}
		
		
/*	Create and store the new sequences.*/
		for(i = 0; i < MAX_PEPTIDE_LENGTH; i++)
		{
			sequence[i] = 0;	/*initialize*/
		}
		sequence[0] = -1;
		cycle = 0;	/*initialize*/
		
		while(Ratchet(aaNum, cycle, sequence, currPtr->peptideLength))
		{
			for(i = 0; i < currPtr->peptideLength; i++)
			{
				peptide[i] = gGapList[sequenceMatrix[i][sequence[i]]];
				aaCorrection[i] = gNodeCorrection[sequenceMatrix[i][sequence[i]]];
				
			}
			peptideLength = currPtr->peptideLength;
			score = currPtr->score;
			nodeValue = currPtr->nodeValue;
			nodeCorrection = currPtr->nodeCorrection;
			gapNum = currPtr->gapNum;
			
			test = GoodSequence(firstSequencePtr, peptide, aaCorrection, peptideLength);
			
			if(test)	/*if different from the original, then store it as a new addition*/
			{
					/*debug*/
				/*	p = 0;
					testPtr = firstSequencePtr;
					while(testPtr != NULL)
					{
						p++;
						testPtr = testPtr->next;
					}*/
				newPtr = LoadSequenceStruct(peptide, peptideLength, score, nodeValue, gapNum, nodeCorrection);
					/*debug*/
				/*	p = 0;
					testPtr = firstSequencePtr;
					while(testPtr != NULL)
					{
						p++;
						testPtr = testPtr->next;
					}*/
				lastPtr->next = newPtr;
				lastPtr = newPtr;
			}
		}
		currPtr = currPtr->next;
	}
	
/*	
	Count the sequences in the list again.
*/
	countTheSeqs = 0;
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		countTheSeqs += 1;
		currPtr = currPtr->next;
	}

	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Sequences expanded to %5ld for qtof score.\n", countTheSeqs);
	}
		 
	return;
}

/*****************************Ratchet************************************************
*
*	This function is similar to RatchetIt that was used for incorporating Edman data.
*	Its used to alter the values in the array sequence.
*
*/

char Ratchet(INT_4 *aaNum, INT_4 cycle, INT_4 *sequence, INT_4 seqLength)
{
	char test = TRUE;
	
	sequence[cycle]++;
	
	if(sequence[cycle] >= aaNum[cycle] && cycle <= seqLength)
	{
		sequence[cycle] = 0;
		cycle++;
		if(cycle > seqLength)
		{
			test = FALSE;
			return(test);
		}
		test = Ratchet(aaNum, cycle, sequence, seqLength);
	}
	return(test);
}
	
	
	
/*************************AddToGapList****************************************
*
*	If a "dipeptide" from a sequence is actually due to three amino acids or more, then
*	it will not be found in the gGapList.  To make this work, I need to add any
*	of these strange "dipeptides" to gGapList.
*/

void	AddToGapList(struct Sequence *firstSequencePtr)
{
	struct Sequence *currPtr;
	INT_4 i, j;
	char test;
	

	
	/*see if any "amino acids" in the sequence do not match any of the gGapList entries*/
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		for(i = 0; i < currPtr->peptideLength; i++)
		{
			test = TRUE;
			for(j = 0; j <= gGapListIndex; j++)
			{
				/*if(currPtr->peptide[i] <= gGapList[j] + gParam.fragmentErr &&
					currPtr->peptide[i] >= gGapList[j] - gParam.fragmentErr)*/
				if(currPtr->peptide[i] == gGapList[j])
				{
					test = FALSE;
					break;
				}
			}
			if(test)	/*if not found, then add this to the list*/
			{
				if(gGapListIndex < MAX_GAPLIST - 1)
				{
					gGapListIndex++;
					gGapList[gGapListIndex] = currPtr->peptide[i];
					gNodeCorrection[gGapListIndex] = 0;
				}
			}
		}
		currPtr = currPtr->next;
	}
	return;
}
	
	
	
/*****************************MakeNewgGapList****************************************
*
*	This function makes a new gGapList where Q and K are differentiated, and the position for
*	'I' becomes 'm' and represents oxidized Met.  Dipeptides are all listed unless they are exactly
*	the same mass (isomeric).  A corresponding gNodeCorrection is used for all gGapList members
*	which represents the value that is added to gGapList to obtain the value of the next decimal
*	place.
*/

void MakeNewgGapList()
{
	INT_4 i, j, k;
	INT_4 lysGlnDiff = (gMonoMass_x100[K] - gMonoMass_x100[Q]) * 1.5;
	INT_4 metPheDiff;
	INT_4 sum, massToAdd, absentFlag;
	REAL_4 correction;
	char delAmAcid, duplicateFlag;
	
	gMonoMass[I] 		= 147.0354;	/*Reset the 'I' position to represent oxidized Met*/
	gAvMass[I] 			= 147.193;
	gNomMass[I] 		= 147;
	gSingAA[I] 			= 'm';
	gMonoMass_x100[I] 	= (gMonoMass[I] * gMultiplier) + 0.5;
	
	metPheDiff 			= (gMonoMass_x100[F] - gMonoMass_x100[I]) * 1.5;
	
	for(i = 0; i < MAX_GAPLIST; i++)
	{
		gGapList[i] = 0;
	}
	gGapListIndex = -1;

	for(i = 0; i < gAminoAcidNumber; i++)	/*Copy the single aa extension masses.*/
	{	
		absentFlag = FALSE;
		if(gParam.aaAbsent[0] != '*') /* Check to see if the AA is on the absent list. */
		{							  /* (We won't add it to the gap list if it is.)   */
			delAmAcid = gParam.aaAbsent[0];	
			j = 0;
			while(delAmAcid != 0 && (delAmAcid >= 'A' && delAmAcid <= 'Y'))
			{
				if(gSingAA[i] == delAmAcid)
				{
					absentFlag = TRUE;
					break;
				}
				j++;
				delAmAcid = gParam.aaAbsent[j];
			}
		}
	
		if(absentFlag || (i == I && gParam.qtofErr >= metPheDiff) || 
			(i == Q && gParam.qtofErr >= lysGlnDiff)) /* Ile and Gln, which are represented by Leu and Lys.*/
		{
			massToAdd = 0;
			correction = 0;
		}
		else if(i == C && gParam.cysMW != 0)
		{
			massToAdd = (gParam.cysMW) + 0.5;	/*Change the mass for cysteine (in case its alkylated).*/
			correction = (gMonoMass[C] * gMultiplier * 10) - (gMonoMass_x100[C] * 10);
		}
		else
		{
			massToAdd = gMonoMass_x100[i];
			correction = (gMonoMass[i] * gMultiplier * 10) - (gMonoMass_x100[i] * 10);
		}
		
		gGapList[i] = massToAdd;
		if(correction >= 0)
		{
			gNodeCorrection[i] = correction + 0.5;
		}
		else
		{
			gNodeCorrection[i] = correction - 0.5;
		}
	}	

	gGapListIndex = gAminoAcidNumber - 1;
	for(i = 0; i < gAminoAcidNumber; i++)	/*Fill in the masses of the 2 AA extensions.*/
	{
		for(j = i; j < gAminoAcidNumber; j++)
		{
			if(gGapList[i] == 0 || gGapList[j] == 0) continue;
			sum = gGapList[i] + gGapList[j];
			/*sum = ((gMonoMass[i] + gMonoMass[j]) * gMultiplier) + 0.5;*/

			correction = ((gMonoMass[i] + gMonoMass[j]) * gMultiplier * 10) - 
							((gMonoMass_x100[i] + gMonoMass_x100[j]) * 10);
			/*correction = (((gMonoMass[i] + gMonoMass[j]) * gMultiplier) * 10) -
							(sum * 10);*/

			duplicateFlag = FALSE;
			for(k = 0; k <= gGapListIndex; k++)
			{
				if(gGapList[k] == sum && (gNodeCorrection[k] <= correction + 1 &&
											gNodeCorrection[k] >= correction - 1))

				{

					/* We already have this mass so don't add it to the list. */
					duplicateFlag = TRUE;
					break;
				}
			}
			
			if(!duplicateFlag)
			{
				gGapListIndex = gGapListIndex + 1;
				gGapList[gGapListIndex] = sum;
				if(correction >= 0)
				{
					gNodeCorrection[gGapListIndex] = correction + 0.5;
				}
				else
				{
					gNodeCorrection[gGapListIndex] = correction - 0.5;
				}
			}
		}
	}

	
	for(i = 0; i <= gGapListIndex; i++)
	{
		if(gNodeCorrection[i] >= 10)
		{
			gNodeCorrection[i] = gNodeCorrection[i] - 10;
			gGapList[i] = gGapList[i] + 1;
		}
		else if(gNodeCorrection[i] <= -10)
		{
			gNodeCorrection[i] = gNodeCorrection[i] + 10;
			gGapList[i] = gGapList[i] - 1;
		}
	}
	
	
	
	return;
}
/*****************************SequenceLengthCalc*************************************
*
*	Calculate the actual sequence length by assuming that any gaps are due to two
*	amino acids rather than counting them as one.
*/

INT_4	SequenceLengthCalc(INT_4 *sequence, INT_4 seqLength)
{
	INT_4 realSeqLength, i, j;
	char test;
	realSeqLength = 0;
	for(i = 0; i < seqLength; i++)
	{
		test = FALSE;	/*test = FALSE when its not a single amino acid, or a two aa extension w/ Pro*/
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] == gGapList[j])
			{
				test = TRUE;
				break;
			}
			if(sequence[i] <= gProPlus[j] + gToleranceNarrow &&
				sequence[i] >= gProPlus[j] - gToleranceNarrow)
			{
				test = TRUE;
				break;
			}
		}
		if(i == 0)	/*don't penalize for a gap at the n-terminal position*/
		{
			test = TRUE;
		}
		if(test)
		{
			realSeqLength++;
		}
		else
		{
			realSeqLength += 2;
		}
	}


	return(realSeqLength);
}

/*****************************SequenceLengthCalc*************************************
*
*	Calculate the actual sequence length by assuming that any gaps are due to two
*	amino acids rather than counting them as one.
*/

INT_4	SequenceLengthCalcNoFudge(INT_4 *sequence, INT_4 seqLength)
{
	INT_4 realSeqLength, residueNumGuess, i, j;
	REAL_4 averageResidueMass = AV_RESIDUE_MASS * gMultiplier;
	char test;
	realSeqLength = 0;
	for(i = 0; i < seqLength; i++)
	{
		test = FALSE;	/*test = FALSE when its not a single amino acid, or a two aa extension w/ Pro*/
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] == gGapList[j])
			{
				test = TRUE;
				break;
			}
		}
		if(test)
		{
			realSeqLength++;
		}
		else
		{
			for(j = gAminoAcidNumber; j < gGapListIndex; j++)
			{
				if(sequence[i] == gGapList[j])
				{
					test = TRUE;
					break;
				}
			}
			if(test)
			{
				realSeqLength += 2;
			}
			else
			{
				residueNumGuess = ((REAL_4)sequence[i] / averageResidueMass) + 0.5;
				realSeqLength = realSeqLength + residueNumGuess;
			}
		}
	}


	return(realSeqLength);
}
/***************************** BoostTheCTerminals ************************************
*
*	For tryptic peptides, the +1 y ions for C-terminal Lys or Arg are boosted for QTof
*	and quadrupole data for +1 or +2 precursors.  For LCQ data, the high mass b ions for
*	cleavage of Lys or Arg are also boosted.  The same boosting is done for C-terminal
*	Lys for Lys-C peptides, Asp or Glu for Glu-C peptides, and Asp for N-terminal Asp-N 
*	peptides.
*/
void BoostTheCTerminals(struct MSData *firstMassPtr)
{
	struct MSData *currPtr;
	INT_4 yArg, yLys, bArg, bLys, yGlu, yAsp, bGlu, bAsp;
	INT_4 bCalStart, yCalStart;
	INT_4 ionNum = 0;
	REAL_4 cTerminalBoost = 2.5;	/*Boost ion intensity by this amount.*/
	REAL_4 averageIntensity = 0;
	REAL_4 lcqBIonIntensity = 0;
	REAL_4 mToAFactor, argInt, lysInt;
	
	if(gParam.chargeState > 2)
		return;	/*don't do anything for precursor charges greater than 2*/
		
/*	Find the average intensity of ions in the list, and calculate the lcqBIonIntensity
	value, which will be twice the average.  This will be the new intensity for high mass
	b ions for loss of K or R in lcq data.*/
	
	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		averageIntensity += currPtr->intensity;
		ionNum++;
		currPtr = currPtr->next;
	}
	if(ionNum == 0)
	{
		printf("BoostTheCTerminals: ionNum = 0");
		exit(1);
	}
	averageIntensity = averageIntensity / ionNum;
	lcqBIonIntensity = averageIntensity * 2;
	
/*	
	Initialize the starting b ion mass.
*/
	bCalStart = gParam.peptideMW - (gParam.modifiedCTerm + 0.5);
	
	/*Correct for average masses.*/
	if(bCalStart > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(bCalStart >= (gParam.monoToAv - gAvMonoTransition))
		{
			mToAFactor = (gParam.monoToAv - bCalStart) / gAvMonoTransition;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(bCalStart >= (gParam.monoToAv - gAvMonoTransition))
	{
		bCalStart = bCalStart * mToAFactor;
	}
	
/*	
	Initialize the starting mass for y ions.
*/
	yCalStart = gParam.modifiedCTerm + (2 * gElementMass_x100[HYDROGEN]) + 0.5;	

	/*Correct for average masses.*/
	if(yCalStart > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(yCalStart >= (gParam.monoToAv - gAvMonoTransition))
		{
			mToAFactor = (gParam.monoToAv - yCalStart) / gAvMonoTransition;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(yCalStart >= (gParam.monoToAv - gAvMonoTransition))
	{
		yCalStart = yCalStart * mToAFactor;
	}
	
/*	Boost ions for tryptic peptides.*/

	if(gParam.proteolysis == 'T')
	{

		yArg = yCalStart + gMonoMass_x100[R];
		yLys = yCalStart + gMonoMass_x100[K];
		bArg = bCalStart - gMonoMass_x100[R];
		bLys = bCalStart - gMonoMass_x100[K];
		
		/*Look to see if both y1 ions are present and pick the most abundant one to boost*/
		argInt = 0;
		lysInt = 0;
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ <= yArg + gParam.fragmentErr 
				&& currPtr->mOverZ >= yArg - gParam.fragmentErr)
			{
				argInt = currPtr->intensity;
			}
			if(currPtr->mOverZ <= yLys + gParam.fragmentErr 
				&& currPtr->mOverZ >= yLys - gParam.fragmentErr)
			{
				lysInt = currPtr->intensity;
			}
			currPtr = currPtr->next;
		}
		
		if(lysInt != 0 && argInt != 0)
		{
			if(lysInt > 2 * argInt)
			{
				argInt = 0;
			}
			if(argInt > 2 * lysInt)
			{
				lysInt = 0;
			}
		}
		
		/*Arg y ion*/
		if(argInt != 0)
		{
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{	
				if(currPtr->mOverZ <= yArg + gParam.fragmentErr 
					&& currPtr->mOverZ >= yArg - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}
				}
				currPtr = currPtr->next;
			}
		}
	
		/*Lys y ion*/
		if(lysInt != 0)
		{
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= yLys + gParam.fragmentErr 
					&& currPtr->mOverZ >= yLys - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}
				}
				currPtr = currPtr->next;
			}
		}
		
		if(gParam.fragmentPattern == 'L')
		{
			/*Look to see if both b ions are present and pick the most abundant one to boost*/
			argInt = 0;
			lysInt = 0;
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= bArg + gParam.fragmentErr 
					&& currPtr->mOverZ >= bArg - gParam.fragmentErr)
				{
					argInt = currPtr->intensity;
				}	
				if(currPtr->mOverZ <= bLys + gParam.fragmentErr 
					&& currPtr->mOverZ >= bLys - gParam.fragmentErr)
				{
					lysInt = currPtr->intensity;
				}
				currPtr = currPtr->next;
			}
		
			if(lysInt != 0 && argInt != 0)
			{
				if(lysInt > 2 * argInt)
				{
					argInt = 0;
				}
				if(argInt > 2 * lysInt)
				{
					lysInt = 0;
				}
			}
			/*Arg b ion*/
			if(argInt != 0)
			{
				currPtr = firstMassPtr;
				while(currPtr != NULL)
				{
					if(currPtr->mOverZ <= bArg + gParam.fragmentErr 
						&& currPtr->mOverZ >= bArg - gParam.fragmentErr)
					{
						currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
						if(currPtr->intensity < averageIntensity)
						{
							currPtr->intensity = averageIntensity;
						}
					}
					currPtr = currPtr->next;
				}
			}
		
			/*Lys b ion*/
			if(lysInt != 0)
			{
				currPtr = firstMassPtr;
				while(currPtr != NULL)
				{
					if(currPtr->mOverZ <= bLys + gParam.fragmentErr 
						&& currPtr->mOverZ >= bLys - gParam.fragmentErr)
					{
						currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
						if(currPtr->intensity < averageIntensity)
						{
							currPtr->intensity = averageIntensity;
						}
					}
					currPtr = currPtr->next;
				}
			}
		}
	}
	
/*	Boost ions for Lys-C peptides.*/

	if(gParam.proteolysis == 'K')
	{

		yLys = yCalStart + gMonoMass_x100[K];
		bLys = bCalStart - gMonoMass_x100[K];
	
		/*Lys y ion*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ <= yLys + gParam.fragmentErr 
				&& currPtr->mOverZ >= yLys - gParam.fragmentErr)
			{
				currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
				if(currPtr->intensity < averageIntensity)
				{
					currPtr->intensity = averageIntensity;
				}
			}
			currPtr = currPtr->next;
		}
		
		if(gParam.fragmentPattern == 'L')
		{
			/*Lys b ion*/
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= bLys + gParam.fragmentErr 
					&& currPtr->mOverZ >= bLys - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}
				}
				currPtr = currPtr->next;
			}
		}
	}
/*	Boost ions for Glu-C peptides.*/

	if(gParam.proteolysis == 'E')
	{

		yGlu = yCalStart + gMonoMass_x100[E];
		yAsp = yCalStart + gMonoMass_x100[D];
		bGlu = bCalStart - gMonoMass_x100[E];
		bAsp = bCalStart - gMonoMass_x100[D];
		
		/*Glu y ion*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{	
			if(currPtr->mOverZ <= yGlu + gParam.fragmentErr 
				&& currPtr->mOverZ >= yGlu - gParam.fragmentErr)
			{
				currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
				if(currPtr->intensity < averageIntensity)
				{
					currPtr->intensity = averageIntensity;
				}
			}
			currPtr = currPtr->next;
		}
	
		/*Asp y ion*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ <= yAsp + gParam.fragmentErr 
				&& currPtr->mOverZ >= yAsp - gParam.fragmentErr)
			{
				currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
				if(currPtr->intensity < averageIntensity)
				{
					currPtr->intensity = averageIntensity;
				}
			}
			currPtr = currPtr->next;
		}
		
		if(gParam.fragmentPattern == 'L')
		{
			/*Glu b ion*/
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= bGlu + gParam.fragmentErr 
					&& currPtr->mOverZ >= bGlu - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}
				}
				currPtr = currPtr->next;
			}
		
			/*Asp b ion*/
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= bAsp + gParam.fragmentErr 
					&& currPtr->mOverZ >= bAsp - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}
				}
				currPtr = currPtr->next;
			}
		}
	}
	
/*	Boost ions for Asp-N peptides.*/

	if(gParam.proteolysis == 'D')
	{
		/*Initialize the starting b ion mass (acetylated, etc).*/
		
		bCalStart = gParam.modifiedNTerm + 0.5;
	
		/*Correct for average masses.*/
		if(bCalStart > gParam.monoToAv)
		{
			mToAFactor = 0;
		}
		else
		{
			if(bCalStart >= (gParam.monoToAv - gAvMonoTransition))
			{
				mToAFactor = (gParam.monoToAv - bCalStart) / gAvMonoTransition;
			}
			else
			{
				mToAFactor = 1;
			}
		}
		mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
		if(bCalStart >= (gParam.monoToAv - gAvMonoTransition))
		{
			bCalStart = bCalStart * mToAFactor;
		}
	
		/*Initialize the starting y ion mass.*/
		yCalStart = gParam.peptideMW + (2 * gElementMass_x100[HYDROGEN]) - (gParam.modifiedNTerm + 0.5);
		
		/*Correct for average masses.*/
		if(yCalStart > gParam.monoToAv)
		{
			mToAFactor = 0;
		}
		else
		{
			if(yCalStart >= (gParam.monoToAv - gAvMonoTransition))
			{
				mToAFactor = (gParam.monoToAv - yCalStart) / gAvMonoTransition;
			}
			else
			{
				mToAFactor = 1;
			}
		}
		mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
		if(yCalStart >= (gParam.monoToAv - gAvMonoTransition))
		{
			yCalStart = yCalStart * mToAFactor;
		}
		
		yAsp = yCalStart - gMonoMass_x100[D];
		bAsp = bCalStart + gMonoMass_x100[D];
			
		/*Asp y ion*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ <= yAsp + gParam.fragmentErr 
				&& currPtr->mOverZ >= yAsp - gParam.fragmentErr)
			{
				currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
				if(currPtr->intensity < averageIntensity)
				{
					currPtr->intensity = averageIntensity;
				}
			}
			currPtr = currPtr->next;
		}
		
		if(gParam.fragmentPattern == 'L')
		{
			/*Asp b ion*/
			currPtr = firstMassPtr;
			while(currPtr != NULL)
			{
				if(currPtr->mOverZ <= bAsp + gParam.fragmentErr 
					&& currPtr->mOverZ >= bAsp - gParam.fragmentErr)
				{
					currPtr->intensity = (currPtr->intensity) * cTerminalBoost;
					if(currPtr->intensity < averageIntensity)
					{
						currPtr->intensity = averageIntensity;
					}	
				}
				currPtr = currPtr->next;
			}
		}
	}
	return;
}
/*********************** StandardDeviationOfTheBYErrors ******************************
*
*	For each peptide, the average error is determined and the standard deviation
*	around this average is determined.  Smaller standard deviations of the error
*	indicate tighter agreement with the proposed sequence.
*/
REAL_4 StandardDeviationOfTheBYErrors(REAL_8 *byError, INT_4 fragNum)
{
	REAL_8 averageError, diffSquared, sumOfDiffSquared, stDevErr;
	REAL_4 stDev;
	INT_4 errorNum, i;
	
	averageError = 0;
	errorNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(byError[i] < 100)
		{
			averageError += byError[i];
			errorNum++;
		}
	}
	
	if(errorNum < 2)
	{
		return(0);
	}
	
	averageError = averageError / errorNum;
	
	sumOfDiffSquared = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(byError[i] < 100)
		{
			diffSquared = byError[i] - averageError;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
	}
	if(errorNum <= 1)
		return(0);
	stDevErr = sumOfDiffSquared / (errorNum - 1);
	stDevErr = sqrt(stDevErr);
	stDev = stDevErr;	/*convert to a REAL_4 to return*/
	return(stDev);
}

/*************************** AssignError *********************************************
*
*	Assigns plus/negative error to ions found as b or y ions.  If an ion has been found
*	as both b and y, then the lower error is saved.
*/
REAL_4 AssignError(REAL_4 currentError, INT_4 calculatedMass, INT_4 observedMass)
{
	REAL_4 newError;
	INT_4 oldMassDiff, newMassDiff;
	
	newMassDiff = abs(calculatedMass - observedMass);
	oldMassDiff = currentError;
	oldMassDiff = abs(oldMassDiff);
	
	newError = observedMass - calculatedMass;
	
	if(currentError != 100)
	{
		if(oldMassDiff < newMassDiff)
		{
			newError = currentError;
		}
	}

	return(newError);
}

/******************************* SingleAA *********************************************
*
*	Returns a FALSE if the amino acid is not single amino acid extension.
*/

char SingleAA(INT_4 aminoAcidMass)
{
	char test = FALSE;
	INT_4 i;
	INT_4 error = 0.4 * gMultiplier;
	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(aminoAcidMass <= gMonoMass_x100[i] + error && 
			aminoAcidMass >= gMonoMass_x100[i] - error)
		{
			test = TRUE;
			break;
		}
	}

	return(test);
}

/************************** RemoveRedundantSequences ***********************************
*
*
*/

struct Sequence *RemoveRedundantSequences(struct Sequence *firstSequencePtr)
{
	struct Sequence *currPtr, *checkThisPtr, *freeMePtr, *previousPtr;
	
	INT_4 i, j, testMass, countTheSeqs;
	INT_4 error = 0.4 * gMultiplier;	/*if error is big, then many unrelated seqs are removed*/
	char test;

	
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		if(currPtr->score != 0)
		{
			checkThisPtr = firstSequencePtr->next;
			while(checkThisPtr != NULL)
			{
				if(checkThisPtr->peptideLength <= currPtr->peptideLength && 
					checkThisPtr->score != 0 && checkThisPtr != currPtr)
				{
					test = TRUE;	/*is TRUE if the two sequences are similar*/
					j = 0;
					for(i = 0; i < checkThisPtr->peptideLength; i++)
					{
						if(j < currPtr->peptideLength)
						{
							if(currPtr->peptide[j] <= checkThisPtr->peptide[i] + error
								&& currPtr->peptide[j] >= checkThisPtr->peptide[i] - error)
							{
								j++;
							}
							else
							{
								if(checkThisPtr->peptide[i] < currPtr->peptide[j])
								{
									test = FALSE;
									break;
								}
								testMass = currPtr->peptide[j];
								test = SingleAA(currPtr->peptide[j]);
								while(j < currPtr->peptideLength && 
									testMass < checkThisPtr->peptide[i] - error)
								{
									j++;
									testMass += currPtr->peptide[j];
									test = SingleAA(currPtr->peptide[j]);
									if(testMass <= checkThisPtr->peptide[i] + error
										&& testMass >= checkThisPtr->peptide[i] - error)
									{
										/*make sure this is not a K = AG situation*/
										if((testMass <= gMonoMass_x100[K] + error 
											&& testMass >= gMonoMass_x100[K] - error) ||
											
											(testMass <= gMonoMass_x100[Q] + error 
											&& testMass >= gMonoMass_x100[Q] - error) ||
											
											(testMass <= gMonoMass_x100[R] + error 
											&& testMass >= gMonoMass_x100[R] - error) ||
											
											(testMass <= gMonoMass_x100[W] + error 
											&& testMass >= gMonoMass_x100[W] - error) ||
											
											(testMass <= gMonoMass_x100[N] + error 
											&& testMass >= gMonoMass_x100[N] - error))
										{
											test = FALSE;
										}
										j++;
										break;
									}
									if(testMass > checkThisPtr->peptide[i] + error)
									{
										test = FALSE;
										break;
									}
								}
								if(test == FALSE)
								{
									break;
								}
							}
						}
					}
					if(test)
					{
						if(checkThisPtr->gapNum == -100)	/*this is signal that sequence 
														was from database*/
						{
							currPtr->score = 0;
						}
						else
						{
							checkThisPtr->score = 0;
						}
					}
				}
				checkThisPtr = checkThisPtr->next;
			}
		
		}
	
		currPtr = currPtr->next;
	}
	
	
	
	/*Free the pointers with scores of zero*/
	while(firstSequencePtr != NULL)
	{
		if(firstSequencePtr->score == 0)
		{
			freeMePtr = firstSequencePtr;
			firstSequencePtr = firstSequencePtr->next;
			free(freeMePtr);
		}
		else
		{
			break;
		}
	}
	
	if(firstSequencePtr != NULL)
	{
		previousPtr = firstSequencePtr;
		currPtr = firstSequencePtr->next;
		while(currPtr != NULL)
		{
			if(currPtr->score == 0)
			{
				freeMePtr = currPtr;
				previousPtr->next = currPtr->next;
				currPtr = currPtr->next;
				free(freeMePtr);
			}
			else
			{
				currPtr = currPtr->next;
				previousPtr = previousPtr->next;
			}
		}
	}
	
	countTheSeqs = 0;
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		countTheSeqs++;
		currPtr = currPtr->next;
	}
	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Scoring %4ld remaining after removing redundant sequences.\n", countTheSeqs);
	}
	if(countTheSeqs == 0)
	{
		printf("RemoveRedundantSequences:  countTheSeqs = 0\n");
		exit(1);
	}

	return(firstSequencePtr);
}
/***************************RevertBackToReals*****************************************
*
*	Divide all of the mass-related gParams values by gMultiplier.  Revert the peak masses
*	back to floats.
*
*/

void RevertBackToReals(struct MSData *firstMassPtr, struct SequenceScore *firstScorePtr)
{
	struct MSData *currPtr;
	struct SequenceScore *currSeqPtr;
	INT_4 i, j;
	
/*	Assign character sequence to peptideSequence field 
(so as to not lose the K/Q differentiation.*/
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		i = 0;
		while(currSeqPtr->peptide[i] != 0)
		{
			for(j = 0; j <= /*gAminoAcidNumber*/ gGapListIndex; j++)
			{
				if(currSeqPtr->peptide[i] == gGapList[j])
				{
					currSeqPtr->peptideSequence[i] = j;
				}
			}
			i++;
		}
		currSeqPtr = currSeqPtr->next;
	}
	
/*	Convert the standard deviations of the b and y errors back.*/
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		currSeqPtr->stDevErr = (currSeqPtr->stDevErr) / gMultiplier;
		currSeqPtr = currSeqPtr->next;
	}

/*	Convert the gParam values back.*/

	gParam.peptideMW   = gParam.peptideMW   / gMultiplier;
	gParam.monoToAv    = gParam.monoToAv    / gMultiplier;
	gParam.peptideErr  = gParam.peptideErr  / gMultiplier;
	gParam.fragmentErr = gParam.fragmentErr / gMultiplier;
	gParam.ionOffset   = gParam.ionOffset   / gMultiplier;
	gParam.cysMW       = gParam.cysMW       / gMultiplier;
	gParam.tagNMass    = gParam.tagNMass    / gMultiplier;
	gParam.tagCMass    = gParam.tagCMass    / gMultiplier;
	gParam.peakWidth   = gParam.peakWidth   / gMultiplier;
	gParam.qtofErr 	   = gParam.qtofErr     / gMultiplier;
	gParam.modifiedNTerm = gParam.modifiedNTerm / gMultiplier;
	gParam.modifiedCTerm = gParam.modifiedCTerm / gMultiplier;
/*	gToleranceWide     = (REAL_4)gToleranceWide     / gMultiplier;
	gToleranceNarrow   = (REAL_4)gToleranceNarrow   / gMultiplier;*/
	
/*	Assign 0 values for gNomMass in places where gGapList is zero (to eliminate using amino
	acids that are absent, or redundant.  I need to do this in order to calculate b and y ions
	in xcorr properly*/
	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		gNomMass[i] = (REAL_4)gGapList[i] / gMultiplier + 0.5;
	}
	
/*	Convert gGapList to nominal values*/
	
	/*for(i = 0; i < gGapListIndex; i++)
	{
		gGapList[i] = gGapList[i] / gMultiplier;
	}
	for(i = 0; i < gGapListIndex; i++)
	{
		for(j = 0; j < gGapListIndex; j++)
		{
			if(i != j)
			{
				if(gGapList[i] != 0 && gGapList[j] != 0)
				{
					if(gGapList[i] == gGapList[j])
					{
						gGapList[i] = 0;
					}
				}
			}
		}
	}*/

/*	Convert the peak masses back.*/

	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		currPtr->mOverZ = (currPtr->mOverZ) / gMultiplier;
		currPtr = currPtr->next;
	}
	
/*	Convert the sequence list.*/
	
	currSeqPtr = firstScorePtr;
	while(currSeqPtr != NULL)
	{
		i = 0;
		while(currSeqPtr->peptide[i] != 0)
		{
			currSeqPtr->peptide[i] = (currSeqPtr->peptide[i]) / gMultiplier;
			i++;
		}
		currSeqPtr = currSeqPtr->next;
	}
	
	if(gAmIHere)
	{
		i = 0;
		while ( gRightSequence[i] != 0)
		{
			gRightSequence[i] = gRightSequence[i] / gMultiplier;
			i++;
		}
	}
	
	return;
}

/***************************RevertTheRevertBackToReals*****************************************
*
*	Change everything back the way it was (I know this is stupid...).
*
*/

void RevertTheRevertBackToReals(struct MSData *firstMassPtr)
{
	struct MSData *currPtr;
	INT_4 i;

/*	Convert the gParam values back.*/

	gParam.peptideMW   = gParam.peptideMW   * gMultiplier;
	gParam.monoToAv    = gParam.monoToAv    * gMultiplier;
	gParam.peptideErr  = gParam.peptideErr  * gMultiplier;
	gParam.fragmentErr = gParam.fragmentErr * gMultiplier;
	gParam.ionOffset   = gParam.ionOffset   * gMultiplier;
	gParam.cysMW       = gParam.cysMW       * gMultiplier;
	gParam.tagNMass    = gParam.tagNMass    * gMultiplier;
	gParam.tagCMass    = gParam.tagCMass    * gMultiplier;
	gParam.peakWidth   = gParam.peakWidth   * gMultiplier;
	gParam.qtofErr 	   = gParam.qtofErr     * gMultiplier;
	gParam.modifiedNTerm = gParam.modifiedNTerm * gMultiplier;
	gParam.modifiedCTerm = gParam.modifiedCTerm * gMultiplier;
	
/*	Assign 0 values for gNomMass in places where gGapList is zero (to eliminate using amino
	acids that are absent, or redundant.  I need to do this in order to calculate b and y ions
	in xcorr properly*/
	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(gGapList[i] == 0)
		{
			gNomMass[i] = 0;
		}
		
		/*gNomMass[i] = (REAL_4)gGapList[i] * gMultiplier;	*/
	}
	
	gMonoMass[9] = 113.08407;
	gMonoMass[6] = 128.05858;
	gSingAA[9] = 'I';

/*	Convert the peak masses back.*/

	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		currPtr->mOverZ = (currPtr->mOverZ) * gMultiplier;
		currPtr = currPtr->next;
	}
	
	return;
}


/******************************* AdjustIonIntensity ******************************************
*
*
*/
void AdjustIonIntensity(INT_4 fragNum, INT_4 *fragIntensity)
{
	REAL_8 averageIntensity = 0;
	REAL_8 intensityDiff = 0;
	REAL_8 summedIntensityDiff = 0;
	REAL_8 standardDev = 0;
	INT_4 confidenceLimitHigh, confidenceLimitLow;
	INT_4 i;
	
	for(i = 0; i < fragNum; i++)
	{
		averageIntensity += fragIntensity[i];
	}
	
	if(fragNum <= 1)
	{
		printf("AdjustIonIntensity:  fragNum = 0\n");
		exit(1);
	}
	
	averageIntensity = (REAL_4)averageIntensity / fragNum + 0.5;
	
	/*Calculate the standard deviation*/
	for(i = 0; i < fragNum; i++)
	{
		intensityDiff = averageIntensity - fragIntensity[i];
		intensityDiff = intensityDiff * intensityDiff;
		summedIntensityDiff += intensityDiff;
	}
	summedIntensityDiff = summedIntensityDiff / (fragNum - 1);
	
	standardDev = sqrt(summedIntensityDiff);
	
	confidenceLimitHigh = averageIntensity + (standardDev * 1.64);	/*1.64 defines 90% conf lim*/
	confidenceLimitLow = averageIntensity - (standardDev * 1.64);	/*1.64 defines 90% conf lim*/
	if(confidenceLimitLow < 0)
	{
		confidenceLimitLow = averageIntensity / 5;
	}
	
	/*Adjust the outliers*/
	for(i = 0; i < fragNum; i++)
	{
		if(fragIntensity[i] > confidenceLimitHigh)
		{
			fragIntensity[i] = confidenceLimitHigh;
		}
	/*	if(fragIntensity[i] < confidenceLimitLow)
		{
			fragIntensity[i] = confidenceLimitLow;
		}*/
	}


	return;
}

/**********************ScoreBYIsotopes*******************************************************
*
*	If the peak width is less than 1.4 Da, then its possible that there are some isotope
*	peaks that were not identified, but are in fact c13 peaks of b or y ions.  Whenever
*	there is a value of one in ionFound and a value of zero in the next spot, if this mass
*	difference is one then the ionFound value for the isotop is given a value of 
*	IONFOUND_ISOTOPE.
*
*/
void ScoreBYIsotopes(REAL_4 *ionFound, INT_4 *fragMOverZ, INT_4 fragNum, INT_4 *ionType)
{
	INT_4 i;
	REAL_4 massDiff, c13MinusC12, oldIonFoundValue, newIonFoundValue;
	
	c13MinusC12 = 1.003354 * gMultiplier;
	
	for(i = 0; i < fragNum - 1; i++)
	{
		if(ionFound[i] != 0)
		{
			massDiff = fragMOverZ[i + 1] - fragMOverZ[i];
			if(massDiff <= c13MinusC12 + gToleranceWide && massDiff >= c13MinusC12 - gToleranceWide)
			{
				oldIonFoundValue = ionFound[i + 1];
				newIonFoundValue = ionFound[i] * IONFOUND_ISOTOPE;
				if(newIonFoundValue > oldIonFoundValue)
				{
					ionFound[i + 1] = newIonFoundValue;
					ionType[i+1] = 15;	/*isotope ions*/
				}
			}
		}
	}

	return;
}

/**********************FreeAllSequenceScore*******************************
*
*	Free linked list of SequenceScore structs.
*
*/
void FreeAllSequenceScore(struct SequenceScore *currPtr)
{
	struct SequenceScore *freeMePtr;
	while(currPtr != NULL)
	{
		freeMePtr = currPtr;
		currPtr = currPtr->next;
		free(freeMePtr);
	}
	return;
}


/********************************AlterIonFound************************************************
*
*	This function multiplies the arrays yFound and bFound, and then checks to see if the
*	corresponding values of fragMOverZ have differences corresponding to the sequence.  In
*	cases where there is a series of three or more ions, those values of ionFound are attenuated
*	by a factor of OVER_USED_IONS.  The cleavageSites is reduced by one for every two such ions.
*
*/
INT_4 AlterIonFound(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				     INT_4 *sequence, INT_4 seqLength, REAL_4 *yFound, 
				     REAL_4 *bFound, INT_4 newCleavageSites)
{
	INT_4 i, j;
	INT_4 yCal, ionsInARow, *theIndexValues, yCalStart;
	REAL_8 mToAFactor;
	char test;

/*	Make room for an array.*/
	theIndexValues = (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(theIndexValues == NULL)
	{
		printf("AlterIonFound: Out of memory.");
		exit(1);
	}

/*	Initialize some variables.*/
	ionsInARow = 0;

/*	
	Initialize the starting mass for y ions.
*/
	yCalStart = gParam.modifiedCTerm + (2 * gElementMass_x100[HYDROGEN]) + 0.5;	

/*	
	Multiply the b and y Found arrays; both can only contain 0 or 1, so multiplying the two
	gives 1 if both are one and zero for other cases.
*/
	for(i = 0; i < fragNum; i++)
	{
		yFound[i] = yFound[i] * bFound[i];
	}
	
/*	
	Calculate the y ions and see if they match with the ions from yFound.
*/
	for(i = (seqLength - 1); i >= 0; i--)
	{
		yCal = yCalStart;
		for(j = (seqLength - 1); j >= i; j--)
		{
			yCal += sequence[j];	/*This is the monoisotopic mass.*/
		}
		
/*	Adjust to average mass if necessary.*/
		if(yCal > gParam.monoToAv)
		{
			mToAFactor = 0;
		}
		else
		{
			if(yCal >= (gParam.monoToAv - gAvMonoTransition))
			{
				mToAFactor = (gParam.monoToAv - yCal) / gAvMonoTransition;
			}
			else
			{
				mToAFactor = 1;
			}
		}
		mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
		if(yCal >= (gParam.monoToAv - gAvMonoTransition))
		{
			yCal = yCal * mToAFactor;
		}
		
/*	Set test to TRUE, which becomes FALSE if a fragment ion matching this yCal is found.*/
		test = TRUE;
		for(j = 0; j < fragNum; j++)
		{
			if(yFound[j] == 1)
			{
				if(fragMOverZ[j] <= yCal + gToleranceWide &&
					fragMOverZ[j] >= yCal - gToleranceWide)
				{
					theIndexValues[ionsInARow] = j;
					ionsInARow++;
					test = FALSE;
					break;
				}
			}
		}
/*	
	test is TRUE if the most recently calculated yCal does not match an observed fragment ion.
	If TRUE, then it first checks to see if there has been a run of more than two ions in
	a row.  If so, it attenuates the ionFound values and resets the cleavageSites value.
*/
		if(test)
		{
			if(ionsInARow > 2)
			{
				for(j = 0; j < ionsInARow; j++)
				{
					ionFound[theIndexValues[j]] = ionFound[theIndexValues[j]] * OVER_USED_IONS;
				}
				ionsInARow = ionsInARow / 2;
				newCleavageSites = newCleavageSites - ionsInARow;
			}
			ionsInARow = 0;
		}
	}
	free(theIndexValues);
	return(newCleavageSites);
}

/*****************************CheckItOut******************************************************
*
*	This function is used to determine if the correct sequence is present.  By setting 
*	gAmIHere to be FALSE, this stuff is always skipped.  If 
*	gAmIHere is TRUE, then this function is activated.
*
*/
BOOLEAN CheckItOutSequenceScore(struct SequenceScore *firstSequencePtr)
{
	struct SequenceScore *currPtr, *correctPtr;
	INT_4 totalSubsequences, rank, correctPeptideLength;
	INT_4 i, j = 0;
	INT_4 lowestScore;
	char test;
	
/*	
	Find the length of the correct sequence.
*/
	correctPeptideLength = 0;
	while ( gRightSequence[correctPeptideLength] != 0)
	{
		correctPeptideLength++;
	}
	
/*	
	Find the lowest score here.
*/
	currPtr = firstSequencePtr;	
	while(currPtr->next != NULL)
	{
		currPtr = currPtr->next;
	}
	lowestScore = currPtr->intensityScore;
	
/*
	Find the correct sequence in the linked list.
*/
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		test = TRUE;
		for(i = 0; i < correctPeptideLength; i++)
		{
			/*if(currPtr->peptide[i] < gRightSequence[i] - gParam.fragmentErr ||
				currPtr->peptide[i] > gRightSequence[i] + gParam.fragmentErr)*/
			if(currPtr->peptide[i] != gRightSequence[i])
			{
				test = FALSE;
				break;
			}
		}
		if(test)	/*If this is the correct subsequence, then..*/
		{
			gAmIHere = TRUE;
			totalSubsequences = 0;
			rank = 1;
			correctPtr = currPtr;
			currPtr = firstSequencePtr;
			while(currPtr != NULL)	/*Count the subsequences and determine the rank.*/
			{
				totalSubsequences++;
				if(currPtr->intensityScore > correctPtr->intensityScore)
				{
					rank++;
				}
				currPtr = currPtr->next;
			}
			
			j++;	/*Stop here in the debugger.*/
			return(TRUE);
			break;
		}
		currPtr = currPtr->next;
	}
	if(test != TRUE)
	{
		j++;	/*Stop here in the debugger for when it doesn't match anymore.*/
		return(FALSE);
	}
}


/*****************************CheckItOut******************************************************
*
*	This function is used to determine if the correct sequence is present.  By setting 
*	gAmIHere to be FALSE, this stuff is always skipped.  If 
*	gAmIHere is TRUE, then this function is activated.
*
*/
BOOLEAN CheckItOut(struct Sequence *firstSequencePtr)
{
	struct Sequence *currPtr, *correctPtr;
	INT_4 totalSubsequences, rank, correctPeptideLength;
	INT_4 i, j = 0;
	INT_4 lowestScore;
	INT_4 error = 0.4 * gMultiplier;
	char test;
	
	if(firstSequencePtr == NULL)
	{
		printf("There were no final sequences to check out.");
		exit(1);
	}
	
/*	
	Find the length of the correct sequence.
*/
	correctPeptideLength = 0;
	while ( gRightSequence[correctPeptideLength] != 0)
	{
		correctPeptideLength++;
	}
	
/*	
	Find the lowest score here.
*/
	currPtr = firstSequencePtr;	
	while(currPtr->next != NULL)
	{
		currPtr = currPtr->next;
	}
	lowestScore = currPtr->score;
	
/*
	Find the correct sequence in the linked list.
*/
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		test = TRUE;
		for(i = 0; i < correctPeptideLength; i++)
		{
			if(currPtr->peptide[i] < gRightSequence[i] - error ||
				currPtr->peptide[i] > gRightSequence[i] + error)
			{
				test = FALSE;
				break;
			}
		}
		if(test)	/*If this is the correct subsequence, then..*/
		{
			gAmIHere = TRUE;
			totalSubsequences = 0;
			rank = 1;
			correctPtr = currPtr;
			currPtr = firstSequencePtr;
			while(currPtr != NULL)	/*Count the subsequences and determine the rank.*/
			{
				totalSubsequences++;
				if(currPtr->score > correctPtr->score)
				{
					rank++;
				}
				currPtr = currPtr->next;
			}
			
			j++;	/*Stop here in the debugger.*/
			return(TRUE);
			break;
		}
		currPtr = currPtr->next;
	}
	if(test != TRUE)
	{
		j++;	/*Stop here in the debugger for when it doesn't match anymore.*/
		return(FALSE);
	}
}



/******************************BCalculator********************************************
*
*	This function calculates singly charged b ion masses. It applies the appropriate
*	average mass correction factor, and returns a INT_4.
*
*/
INT_4 BCalculator(INT_4 i, INT_4 *sequence, INT_4 bCalStart, INT_4 bCalCorrection)
{
	INT_4 bCal = bCalStart;
	INT_4 j, corrections[10], residueCount, k, maxCorrection, minCorrection, correctionSpread, avCorrection;
	INT_4 nodeCorrection = bCalCorrection;
	REAL_8 mToAFactor;
	char test;
	
	for(j = 0; j < 10; j++)
	{
		corrections[j] = 0;
	}

	for(j = 0 ; j < i; j++)
	{
		bCal += sequence[j];
	}
	
	for(j = 0; j < i; j++)	/*determine the correction*/
	{
		test = TRUE;
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(sequence[j] == gGapList[k])
			{
				nodeCorrection += gNodeCorrection[k];
				test = FALSE;
			}
		}
		if(test)	/*test is true if this is a dipeptide*/
		{
			residueCount = 0;
			for(k = gAminoAcidNumber; k <= gGapListIndex; k++)
			{
				if(sequence[j] == gGapList[k])
				{
					corrections[residueCount] = gNodeCorrection[k];
					residueCount++;
					if(residueCount >= 10)
					{
						printf("LutefiskScore: residueCount exceeds 10.\n");
						exit(1);
					}
				}
			}
			if(residueCount <= 1)	/*if not in gGapList, then no correction is added, cuz
									corrections array initialized to zero.*/
			{
				nodeCorrection += corrections[0];
			}
			else
			{
				maxCorrection = corrections[0];
				minCorrection = corrections[0];
				for(k = 1; k < residueCount; k++)
				{
					if(corrections[k] > maxCorrection)
					{
						maxCorrection = corrections[k];
					}
					if(corrections[k] < minCorrection)
					{
						minCorrection = corrections[k];
					}
				}
				correctionSpread = maxCorrection - minCorrection;
				if(correctionSpread <= 5)
				{
					avCorrection = 0;
					for(k = 0; k < residueCount; k++)
					{
						avCorrection += corrections[k];
					}
					if(avCorrection > 0)
					{
						avCorrection = ((REAL_4)avCorrection / residueCount) + 0.5;
					}
					else
					{
						avCorrection = ((REAL_4)avCorrection / residueCount) - 0.5;
					}
					nodeCorrection += avCorrection;
				}
			}
		}
	}
	
	
	
	if(nodeCorrection >= 5)
	{
		nodeCorrection = ((REAL_4)nodeCorrection / 10) + 0.5;
		bCal += nodeCorrection;
	}
	else if(nodeCorrection <= -5)
	{
		nodeCorrection = ((REAL_4)nodeCorrection / 10) - 0.5;
		bCal += nodeCorrection;
	}

	
	
	
	
	
	
	if(bCal > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(bCal >= (gParam.monoToAv - gAvMonoTransition))
		{
			mToAFactor = (gParam.monoToAv - bCal) / gAvMonoTransition;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(bCal >= (gParam.monoToAv - gAvMonoTransition))
	{
		bCal = bCal * mToAFactor;
	}

	

	return(bCal);
}

/******************************YCalculator********************************************
*
*	This function calculates singly charged y ion masses.  It applies the appropriate
*	average mass correction factor, and returns a INT_4.
*
*/
INT_4 YCalculator(INT_4 i, INT_4 *sequence, INT_4 seqLength, INT_4 yCalStart, INT_4 yCalCorrection)
{
	INT_4 yCal = yCalStart;
	INT_4 j, k, nodeCorrection = yCalCorrection;
	INT_4 residueCount, correctionSpread, corrections[10], maxCorrection, minCorrection;
	INT_4 avCorrection;
	REAL_8 mToAFactor;
	char test;
	
	for(j = 0; j < 10; j++)
	{
		corrections[j] = 0;
	}
	for(j = i; j < seqLength; j++)
	{
		yCal += sequence[j];
	}
	for(j = i; j < seqLength; j++)
	{
		test = TRUE;
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(sequence[j] == gGapList[k])
			{
				nodeCorrection += gNodeCorrection[k];
				test = FALSE;
			}
		}
		if(test)	/*test is true if this is a dipeptide*/
		{
			residueCount = 0;
			for(k = gAminoAcidNumber; k <= gGapListIndex; k++)
			{
				if(sequence[j] == gGapList[k])
				{
					corrections[residueCount] = gNodeCorrection[k];
					residueCount++;
					if(residueCount >= 10)
					{
						printf("LutefiskScore: residueCount exceeds 10.\n");
						exit(1);
					}
				}
			}
			if(residueCount <= 1)	/*if not in gGapList, then nothing is added to 
									nodecorrection, since corrections array initialized
									to zero*/
			{
				nodeCorrection += corrections[0];
			}
			else
			{
				maxCorrection = corrections[0];
				minCorrection = corrections[0];
				for(k = 1; k < residueCount; k++)
				{
					if(corrections[k] > maxCorrection)
					{
						maxCorrection = corrections[k];
					}
					if(corrections[k] < minCorrection)
					{
						minCorrection = corrections[k];
					}
				}
				correctionSpread = maxCorrection - minCorrection;
				if(correctionSpread <= 5)
				{
					avCorrection = 0;
					for(k = 0; k < residueCount; k++)
					{
						avCorrection += corrections[k];
					}
					if(avCorrection > 0)
					{
						avCorrection = ((REAL_4)avCorrection / residueCount) + 0.5;
					}
					else
					{
						avCorrection = ((REAL_4)avCorrection / residueCount) - 0.5;
					}
					nodeCorrection += avCorrection;
				}
			}
		}
	}
	if(nodeCorrection >= 5)
	{
		nodeCorrection = ((REAL_4)nodeCorrection / 10) + 0.5;
		yCal += nodeCorrection;
	}
	else if(nodeCorrection <= -5)
	{
		nodeCorrection = ((REAL_4)nodeCorrection / 10) - 0.5;
		yCal += nodeCorrection;
	}
		
	
	if(yCal > gParam.monoToAv)
	{
		mToAFactor = 0;
	}
	else
	{
		if(yCal >= (gParam.monoToAv - gAvMonoTransition))
		{
			mToAFactor = (gParam.monoToAv - yCal) / gAvMonoTransition;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	if(yCal >= (gParam.monoToAv - gAvMonoTransition))
	{
		yCal = yCal * mToAFactor;
	}
	
	return(yCal);
}

/***************************TwoAAExtFinder*********************************
*
*	Returns a char value of TRUE or FALSE depending on if the extension contains two
*	amino acids (TRUE) or not (FALSE).
*
*/
char TwoAAExtFinder(INT_4 *sequence, INT_4 i)
{
	INT_4 j;
	char twoAAExtension;
	
	twoAAExtension = TRUE;
	j = 0;
	while(twoAAExtension && j < gAminoAcidNumber)
	{
		if((sequence[i] <= (gMonoMass_x100[j] + gToleranceNarrow)) 
			&& (sequence[i] >= (gMonoMass_x100[j] - gToleranceNarrow)))
		{
			twoAAExtension = FALSE;
		}
		j++;
	}
	return(twoAAExtension);
}

/******************************FindNOxMet****************************************
*
*	This function counts the number of oxidized Mets in a sequence, which is used
*	to keep track of the number of ox met in each b ion.
*
*/
INT_4	FindNOxMet(INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, oxMetCount;
	
	if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0)	/*mass accuracy sufficient
																to determine oxMet*/
	{			
		oxMetCount = 0;
		for(i = 0; i < seqLength; i++)
		{
			if(sequence[i] >= gMonoMass_x100[9] - gToleranceNarrow
			&& sequence[i] <= gMonoMass_x100[9] + gToleranceNarrow)
			{
				oxMetCount++;
			}
		}
	}
	else	/*mass accuracy not sufficient to differentiate oxMet from Phe*/
	{
		oxMetCount = 0;
		for(i = 0; i < seqLength; i++)
		{
			if(sequence[i] >= gMonoMass_x100[F] - gToleranceNarrow
			&& sequence[i] <= gMonoMass_x100[F] + gToleranceNarrow)
			{
				oxMetCount++;
			}
		}
	}
	return(oxMetCount);
}

/******************************FindNCharge***************************************
*
*	This function counts the number of basic residues in a sequence, and returns
*	a INT_4 that contains this number.
*
*/
INT_4 FindNCharge(INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, j, nChargeCount;
	
	nChargeCount = 1;
	for(i = 0; i < seqLength; i++)
	{
		if(((sequence[i] <= (gMonoMass_x100[R] + gToleranceWide)) &&
			(sequence[i] >= (gMonoMass_x100[R] - gToleranceWide))) || 
			(sequence[i] <= (gMonoMass_x100[H] + gToleranceWide) &&
			(sequence[i] >= (gMonoMass_x100[H] - gToleranceWide))) ||
			(sequence[i] <= (gMonoMass_x100[K] + gToleranceWide) &&
			(sequence[i] >= (gMonoMass_x100[K] - gToleranceWide))))
		{
			nChargeCount += 1;
		}
	}
	
	for(i = 0; i < seqLength; i++)	/*Here I look for two amino acid extensions containing
									Arg, His, or Lys.*/
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(((sequence[i]  <=  (gArgPlus[j] + gToleranceWide)) &&
				(sequence[i]  >=  (gArgPlus[j] - gToleranceWide))) ||
				((sequence[i] <=  (gHisPlus[j] + gToleranceWide)) &&
				(sequence[i]  >=  (gHisPlus[j] - gToleranceWide))) ||
				((sequence[i] <=  (gLysPlus[j] + gToleranceWide)) &&
				(sequence[i]  >=  (gLysPlus[j] - gToleranceWide))))
			{
				nChargeCount += 1;
			}
		}
	}
	
	return(nChargeCount);
}


/******************************InitLutefiskScore***************************************
*
*	This function assigns space to some arrays, counts the sequences in the final list of
*	completed sequences, puts the sequence tag back into the list of sequences, determines
*	if there is a C-terminal Lys or Arg for tryptic peptides, and multiplies several mass
*	variables by 100 so that integers can be used rather than REAL_4s.
*
*/
struct Sequence *InitLutefiskScore(INT_4 *sequence, INT_4 *fragMOverZ, INT_4 *fragIntensity,
						 REAL_4 *ionFound, REAL_4 *ionFoundTemplate, INT_4 *countTheSeqs,
						 struct Sequence *firstSequencePtr,
						 REAL_4 *yFound, REAL_4 *bFound, struct MSData *firstMassPtr, 
						 INT_4 *charSequence, REAL_8 *byError, INT_4 *ionType)
{
	struct Sequence *currSeqPtr;
	INT_4 i;
	REAL_4 lowMassIonConversion;
	
/*	
	Check that arrays are ok.
*/
	if(yFound == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(bFound == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(sequence == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(charSequence == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(fragMOverZ == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(fragIntensity == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(ionFound == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(ionType == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(byError == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	if(ionFoundTemplate == NULL)
	{
		printf("InitLutefiskScore:  Out of memory");
		exit(1);
	}
	
/*	
	Count the sequences in the list, print the number to the console, and exit if no sequences.
*/
	*countTheSeqs = 0;
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		*countTheSeqs += 1;
		currSeqPtr = currSeqPtr->next;
	}
	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Scoring %4ld completed sequences.\n", *countTheSeqs);
	}
	if(*countTheSeqs == 0)
	{
		printf("InitLutefiskScore:  countTheSeqs = 0\n");
		exit(1);
	}

/*	
	Insert the sequence tag back into the sequence before scoring it.
*/
	if(gParam.tagSequence[0] != '*')
	{
		firstSequencePtr = AddTagBack(firstSequencePtr);
	}
	
/*	
	Check to see if there is an ion at 147 or 175, indicating a tryptic C-terminus.
	If gTrypticCterm is TRUE then something is done in "MassageScores".
*/
	if(gParam.proteolysis == 'T')
	{
		if(gParam.fragmentPattern == 'Q' || gParam.fragmentPattern == 'T')
		{
			CTerminalLysOrArg(firstMassPtr);
		}
		else
		{
			gTrypticCterm = TRUE;	/*Its difficult to tell if Lys or Arg is present in 
									  LCQ data.*/
		}
	}
	
/*	
	Setup the global arrays gCysPlus, gArgPlus, gHisPlus, and gLysPlus so that they contain
	the appropriate values for the type of cysteine alkyl group.  But don't do this again, cuz 
	the numbers get screwy.
*/
	if(gFirstTimeThru)
	{
		lowMassIonConversion = (REAL_4)gMultiplier / 10000;
		
		for(i = 0; i < gAminoAcidNumber; i++)
		{
			if(gGapList[i] == 0)
			{	
				gArgPlus[i] = 0;
				gHisPlus[i] = 0;
				gLysPlus[i] = 0;
				gCysPlus[i] = 0;
				gGlnPlus[i] = 0;
				gGluPlus[i] = 0;
				gProPlus[i] = 0;
			}
			else
			{
				/*gArgPlus[i] = (REAL_4)gArgPlus[i] * lowMassIonConversion + 0.5;
				gHisPlus[i] = (REAL_4)gHisPlus[i] * lowMassIonConversion + 0.5;
				gLysPlus[i] = (REAL_4)gLysPlus[i] * lowMassIonConversion + 0.5;
				gCysPlus[i] = (REAL_4)gCysPlus[i] * lowMassIonConversion + 0.5;
				gGlnPlus[i] = (REAL_4)gGlnPlus[i] * lowMassIonConversion + 0.5;		
				gGluPlus[i] = (REAL_4)gGluPlus[i] * lowMassIonConversion + 0.5;	
				gProPlus[i] = (REAL_4)gProPlus[i] * lowMassIonConversion + 0.5;*/
					
				gArgPlus[i] = gMonoMass_x100[R] + gMonoMass_x100[i];
				gHisPlus[i] = gMonoMass_x100[H] + gMonoMass_x100[i];
				gLysPlus[i] = gMonoMass_x100[K] + gMonoMass_x100[i];
				gCysPlus[i] = gMonoMass_x100[C] + gMonoMass_x100[i];
				gGlnPlus[i] = gMonoMass_x100[Q] + gMonoMass_x100[i];		
				gGluPlus[i] = gMonoMass_x100[E] + gMonoMass_x100[i];
				gProPlus[i] = gMonoMass_x100[P] + gMonoMass_x100[i];
			}
		}
	
		gArgPlus[C] = gParam.cysMW + gMonoMass_x100[R];
		gHisPlus[C] = gParam.cysMW + gMonoMass_x100[H];
		gLysPlus[C] = gParam.cysMW + gMonoMass_x100[K];
		gGlnPlus[C] = gParam.cysMW + gMonoMass_x100[Q];
		gGluPlus[C] = gParam.cysMW + gMonoMass_x100[E];
		gCysPlus[C] = gParam.cysMW + gParam.cysMW;
		gProPlus[C] = gParam.cysMW + gMonoMass_x100[P];
	
		for(i = 0; i < gAminoAcidNumber; i++)
		{
			if(gGapList[i] == 0)
			{
				gCysPlus[i] = gParam.cysMW + gMonoMass_x100[i];
			}
		}
	}
	
	return(firstSequencePtr);
}

/********************************* HighMOverZFilter ***********************************
*
*
*
*/
void HighMOverZFilter(struct Sequence *firstSequencePtr, INT_4 *fragMOverZ, 
						INT_4 *fragIntensity, INT_4 *countTheSeqs, 
						INT_4 *sequence, INT_4 fragNum)
{
	INT_4 precursor, i, j, seqLength, seqLimit;
	INT_4 *highMZFrags, *highMZInts, highMZNum;
	INT_4 *tempMZFrags, *tempMZInts, tempMZNum, intScoreCutoff;
	INT_4 maxNumOfIons, greatestInt, greatestIntIndex;
	REAL_4 totalIntensity, score, averageIntensity, scoreCutoff;
	REAL_4 obsdIsotopeRatio, calcIsotopeRatio;
	REAL_4 C13minusC12 = 1.00335 * gMultiplier;
	char test;
	
	struct Sequence *currSeqPtr, *previousSeqPtr;
	
	/*Assign some space for these arrays.*/
	highMZFrags = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(highMZFrags == NULL)
	{
		printf("HighMOverZFilter:  Out of memory.");
		exit(1);
	}
	
	highMZInts = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(highMZInts == NULL)
	{
		printf("HighMOverZFilter:  Out of memory.");
		exit(1);
	}
	tempMZFrags = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(highMZFrags == NULL)
	{
		printf("HighMOverZFilter:  Out of memory.");
		exit(1);
	}
	
	tempMZInts = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(highMZInts == NULL)
	{
		printf("HighMOverZFilter:  Out of memory.");
		exit(1);
	}
	/*Initialize variables.*/
	precursor = (gParam.peptideMW + 
				(gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	totalIntensity = 0;
	/*seqLimit is the num of seqs required before the weeding is stopped*/
	if(gParam.chargeState <= 2)
	{
		seqLimit = 50;
	}
	else
	{
		seqLimit = 75;
	}
	/*set a limit on the num of ions*/
	if(gParam.fragmentPattern == 'L')	/*LCQ data has high mass b and y*/
	{
		maxNumOfIons = ((REAL_4)gParam.peptideMW / (2 * gAvResidueMass)) * 1.0 + 0.5; 
	}
	else	/*TSQ and QTOF data generally has only y ions (fewer ions to look at)*/
	{
		maxNumOfIons = ((REAL_4)gParam.peptideMW / (2 * gAvResidueMass)) + 0.5;
	}
	
/*	Generate the list of ions greater than the precursor m/z.*/
	averageIntensity = 0;
	highMZNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(fragMOverZ[i] > precursor + gToleranceWide * 2)
		{
			averageIntensity += fragIntensity[i];
			highMZNum++;
		}
	}
	
	if(highMZNum == 0)
	{
		printf("No high mass ions found for the high m/z filter.\n");
		free(highMZFrags);
		free(highMZInts);
		free(tempMZFrags);
		free(tempMZInts);
		return;
	}
	
	averageIntensity = averageIntensity / highMZNum;
	averageIntensity = averageIntensity / 2;	/*2 is arbitrary; averageIntensity is used below
												in order to select only the more abundant ions;
												if this threshold is high then you might consider
												raising the percentage up from 90% (below)*/
	
	tempMZNum = 0;
	for(i = 0; i < fragNum; i++)
	{
		if(fragMOverZ[i] > precursor + 5 * gMultiplier && fragIntensity[i] > averageIntensity)
		{
			test = TRUE;	/*make sure I don't use an isotope peak*/
			if(tempMZNum != 0)
			{
				for(j = 0; j < tempMZNum; j++)
				{
					if(tempMZFrags[j] + C13minusC12 >= fragMOverZ[i] - gToleranceWide * 2 
						&& tempMZFrags[j] + C13minusC12 <= fragMOverZ[i] + gToleranceWide * 2)
					{
						if(fragIntensity[i] == 0 || tempMZFrags[j] == 0)
						{
							printf("HighMOverZFilter: divide by zero");
							exit(1);
						}
						obsdIsotopeRatio = (REAL_4)tempMZInts[j] / (REAL_4)fragIntensity[i];
						calcIsotopeRatio = (REAL_4)(1800 * gMultiplier) / (REAL_4)tempMZFrags[j];
						if(obsdIsotopeRatio > calcIsotopeRatio - 0.3 ||
							(obsdIsotopeRatio <= 1.2 && obsdIsotopeRatio >= 0.8))
						{
							test = FALSE;
						}
						/*My LCQ gives funny isotope ratios, for some reason*/
						if(gParam.fragmentPattern == 'L')
						{
							if(obsdIsotopeRatio * 2 > calcIsotopeRatio - 0.3)
							{
								test = FALSE;
							}
						}
					}
				}
			}
			if(test)
			{
				tempMZFrags[tempMZNum] = fragMOverZ[i];
				tempMZInts[tempMZNum]  = fragIntensity[i];
				tempMZNum++;
			}
		}
	}
	
	if(tempMZNum <= maxNumOfIons)
	{
		for(i = 0; i < tempMZNum; i++)
		{
			highMZFrags[i] = tempMZFrags[i];
			highMZInts[i] = tempMZInts[i];
		}
		highMZNum = tempMZNum;
	}
	else
	{
		highMZNum = 0;
		for(i = 0; i < maxNumOfIons; i++)
		{
			greatestInt = tempMZInts[0];
			greatestIntIndex = 0;
			for(j = 0; j < tempMZNum; j++)
			{
				if(tempMZInts[j] > greatestInt)
				{
					greatestInt = tempMZInts[j];
					greatestIntIndex = j;
				}
			}
			tempMZInts[greatestIntIndex] = -1 * tempMZInts[greatestIntIndex];
		}
		highMZNum = 0;
		for(i = 0; i < tempMZNum; i++)
		{
			if(tempMZInts[i] < 0)
			{
				highMZFrags[highMZNum] = tempMZFrags[i];
				highMZInts[highMZNum] = -1 * tempMZInts[i];
				highMZNum++;
			}
		}
	}

/*	Calculate the total ion intensity for this high m/z ion list.*/
	for(i = 0; i < highMZNum; i++)
	{
		totalIntensity += highMZInts[i];
	}
	
	if(totalIntensity == 0)
	{
		printf("HighMOverZFilter: totalIntensity = 0\n");
		free(highMZFrags);
		free(highMZInts);
		free(tempMZFrags);
		free(tempMZInts);
		return;
	}
	
/*	Calculate the initial scoreCutoff value by finding the lowest intensity ion, and assume 
	that the best sequences might not account for that one.*/
	
	scoreCutoff = totalIntensity;
	for(i = 0; i < highMZNum; i++)
	{
		if(highMZInts[i] < scoreCutoff)
		{
			scoreCutoff = highMZInts[i];
		}
	}
	scoreCutoff = (totalIntensity - scoreCutoff) / totalIntensity;
	scoreCutoff = scoreCutoff * 100;

/*	Start looking at each sequence.*/
	while(scoreCutoff > 0)
	{
		currSeqPtr = firstSequencePtr;
		while(currSeqPtr != NULL)
		{
			/*if(currSeqPtr->score == 629 && currSeqPtr->nodeValue == 14603)
			{
				i++;
				i++;
			} for debugging*/
			
			/*Load the sequence array.*/
			LoadSequence(sequence, &seqLength, currSeqPtr);
			
			/*	Return a REAL_4 value that corresponds to the intensity only score 
			using b, y and -17/18*/
			score = AssignHighMZScore(highMZNum, highMZFrags, highMZInts, totalIntensity,
										sequence, seqLength);
										
			/*Skip sequences that have Pro in third position, for LCQ data since these tend to
			give lower abudance high mass ions.*/
			
			/*score = ProInThirdPosition(score, sequence, seqLength);*/
			
			/*currSeqPtr->score is a INT_4 and intScore is a REAL_4 that is less than 1,
			so I multiply it by 100 to get an int (%).*/
			if(score >= 0)
			{							
				currSeqPtr->score = score * 100 + 0.5;
			}
			else
			{
				currSeqPtr->score = -1; /*has a proline in the third position*/
			}	
			
			/*Give sequences with scores less than 90% a zero as a signal to trash it later.*/
			intScoreCutoff = scoreCutoff + 0.5;
			if(currSeqPtr->score < intScoreCutoff && currSeqPtr->score > 0)	
			{
				currSeqPtr->score = 0;
			}
	
			currSeqPtr = currSeqPtr->next;
		}
		
		*countTheSeqs = 0;
		currSeqPtr = firstSequencePtr;
		while(currSeqPtr != NULL)
		{
			if(currSeqPtr->score > 0)
			{
				*countTheSeqs += 1;
			}
			currSeqPtr = currSeqPtr->next;
		}
		if(*countTheSeqs < seqLimit)
		{
			scoreCutoff = scoreCutoff - 5;	/*reduce the cutoff by 5%*/
			currSeqPtr = firstSequencePtr;
			while(currSeqPtr != NULL)
			{
				currSeqPtr->score = 1;
				currSeqPtr = currSeqPtr->next;
			}
		}
		else
		{
			if(gCorrectMass && gParam.fMonitor)
			{
				printf("The cutoff for the high m/z ions is %4.1f percent.\n", scoreCutoff);
			}
			scoreCutoff = -1;
		}
	}
/*	Get rid of sequences that have a score of zero, but don't free the firstSequencePtr.*/
	currSeqPtr = firstSequencePtr->next;
	previousSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)	
	{
		if(currSeqPtr->score == 0)
		{
			previousSeqPtr->next = currSeqPtr->next;
			free(currSeqPtr);
			currSeqPtr = previousSeqPtr->next;
		}
		else
		{
			previousSeqPtr = currSeqPtr;
			currSeqPtr = currSeqPtr->next;
		}
	}
	
/*	
	Count the sequences in the list again.
*/
	*countTheSeqs = 0;
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		*countTheSeqs += 1;
		currSeqPtr = currSeqPtr->next;
	}

	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Scoring %5ld sequences after the high m/z filter.\n", *countTheSeqs);
	}
	
	/*Free the arrays.*/
	free(highMZFrags);
	free(highMZInts);
	free(tempMZFrags);
	free(tempMZInts);
	
	return;
}

/************************************TossTheLosers*************************************
*
*	The best tryptic peptide sequences have contiguous series of y and, to a lesser extent,
*	b ions.  Sequences that seem to be derived from wild mixtures of odd ions are weeded out
*	at this point.  This function does not operate on non-tryptic peptides.
*
*/
void TossTheLosers(struct Sequence *firstSequencePtr, REAL_4 *ionFoundTemplate, INT_4 fragNum,
					INT_4 *fragMOverZ, INT_4 *fragIntensity, INT_4 intensityTotal,
					REAL_4 *ionFound, INT_4 *countTheSeqs, INT_4 *sequence)
{
	struct Sequence *currSeqPtr, *previousSeqPtr;
	INT_4 seqLength, i, cleavageSites;
	INT_4 highestBYScore = 0, averageBYScore = 0;
	REAL_4 intScore;
	
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)	/*Assign a BY intensity score to each sequence.*/
	{
		if(currSeqPtr->score == 354 && currSeqPtr->nodeValue == 134773
			&&currSeqPtr->nodeCorrection == -5)
		{
			i++;	
			i++;
		} /*for debugging*/
		
		
		LoadSequence(sequence, &seqLength, currSeqPtr);	/*Also used below.*/
			
		for(i = 0; i < fragNum; i++)	/*Initialize for each sequence.*/
		{
			ionFound[i] = ionFoundTemplate[i];
		}
			
		cleavageSites = FindBYIons(ionFound, fragNum, fragMOverZ, sequence, 
			                       seqLength);	/*Note: This is not FindABYIons!*/
		ProlineInternalFrag(ionFound, fragMOverZ, sequence, seqLength, fragNum);
		intScore = BYIntensityScorer(fragIntensity, ionFound, cleavageSites, fragNum, 
									seqLength, intensityTotal);
		/*save the subsequence score to nodeValue*/
		/*currSeqPtr->nodeValue = currSeqPtr->score;*/
		currSeqPtr->score = (intScore * 1000);	/*currSeqPtr->score is a INT_4 and
													intScore is a REAL_4 that is less than 1,
													so I multiply it by 1000 to get an int.*/
		if(currSeqPtr->score > highestBYScore)
		{
			highestBYScore = currSeqPtr->score;
		}
		averageBYScore += (intScore * 1000);
		currSeqPtr = currSeqPtr->next;
	}

	if(*countTheSeqs == 0)
	{
		printf("TossTheLosers: countTheSeqs = 0\n");
		return;
	}
	averageBYScore = averageBYScore / *countTheSeqs;	/*countTheSeqs is the number of sequences.*/

/*
	Make the average BYScore be an average of the highest and the average score, ie, I only 
	keep the highest quartile.
*/
	averageBYScore = (averageBYScore + highestBYScore) / 2.05;	
																
/*
	For sequences with below average score, I reassign the score field to zero.
*/
	currSeqPtr = firstSequencePtr->next;	
	while(currSeqPtr != NULL)
	{
		if(currSeqPtr->score < averageBYScore)
		{
			currSeqPtr->score = 0;
		}
		currSeqPtr = currSeqPtr->next;
	}

/*
	Now I get rid of the zero score sequences.
*/
	currSeqPtr = firstSequencePtr->next;
	previousSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)	
	{
		if(currSeqPtr->score == 0)
		{
			previousSeqPtr->next = currSeqPtr->next;
			free(currSeqPtr);
			currSeqPtr = previousSeqPtr->next;
		}
		else
		{
			previousSeqPtr = currSeqPtr;
			currSeqPtr = currSeqPtr->next;
		}
	}
	
/*	
	Count the sequences in the list again.
*/
	*countTheSeqs = 0;
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		*countTheSeqs += 1;
		currSeqPtr = currSeqPtr->next;
	}

	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Scoring %4ld sequences following the b and y filter.\n", *countTheSeqs);
	}
	
	return;
	
}


/*******************************CTerminalLysOrArg***************************************
*
*	This function will change the value of gTrypticCTerm (a char) to TRUE if it finds
*	an ion at m/z 147 or 175, and if its been specified that the peptide is derived from
*	tryptic digestion.
*
*/
void CTerminalLysOrArg(struct MSData *firstMassPtr)
{
	REAL_4 lysYIon, argYIon;
	struct MSData *currPtr;
	
	lysYIon = gMonoMass_x100[K] + 3 * gElementMass_x100[HYDROGEN] + gElementMass_x100[OXYGEN];
	argYIon = gMonoMass_x100[R] + 3 * gElementMass_x100[HYDROGEN] + gElementMass_x100[OXYGEN];
	
	currPtr = firstMassPtr;
	gTrypticCterm = FALSE;
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > argYIon + gToleranceWide)
		{
			break;
		}
		if(currPtr->mOverZ <= (lysYIon + gToleranceWide) &&
			currPtr->mOverZ >= (lysYIon - gToleranceWide))
		{
			gTrypticCterm = TRUE;
			break;
		}
		if(currPtr->mOverZ <= (argYIon + gToleranceWide) &&
			currPtr->mOverZ >= (argYIon - gToleranceWide))
		{
			gTrypticCterm = TRUE;
			break;
		}
		currPtr = currPtr->next;
	}
	return;
}

/******************************BYIntensityScorer*********************************
*
*	IntensityScorer inputs fragIntensity (the ion intensities), ionFound (the array that 
*	is indexed the same as fragIntensity and contains 1 for ions that have been identified 
*	and 0 for those that were not), cleavageSites (the number of times either a b or y ion 
*	of any charge was found to delineate the sequence), and fragNum (the number of fragment
*	ions in the CID data).  This function returns a REAL_4 value (ranging from zero to one) 
*	corresponding to the fraction of the ion current that can be identified times a multiplier 
*	that reflects the idea that correct sequences are usually delineated by series of either 
*	b or y ions.
*/
REAL_4 BYIntensityScorer(INT_4 *fragIntensity, REAL_4 *ionFound, INT_4 cleavageSites, 
						INT_4 fragNum, INT_4 seqLength, INT_4 intensityTotal)
{
	REAL_4 intScore = 0;
	REAL_4 numScore = 0;
	REAL_4 attenuation;
	INT_4 i;
	INT_4 intNum = 0;
	
/*	
*	Initialize attenuation, which is a fractional multiplier that reflects the number of 
*	times either b or y ions delineate
*	the proposed sequence.
*/
	attenuation = cleavageSites;
	if(seqLength - 1 == 0)
	{
		printf("BYIntensityScorer: seqLength - 1 = 0\n");
		return(intScore);
	}
	attenuation = attenuation / (seqLength - 1);
	
/*	
	Add up the intensity that has been identified, and count the ions.
*/
	for(i = 0; i < fragNum; i++)
	{
		if(ionFound[i] != 0)
		{
			intScore += fragIntensity[i] * ionFound[i];
			intNum++;
		}
	}
	
/*	
	Divide by the total ion intensity (which does not include the precursor ion region).
*/
	if(intensityTotal == 0)
	{
		printf("BYIntensityScorer:  intensityTotal = 0\n");
		exit(1);
	}
	intScore = intScore / intensityTotal;
	
/*	
	Attenuate the score so that sequences with lots of b and y ions are favored.
*/
	intScore = ((INTENSITY_WEIGHT * intScore) + (ATTENUATION_WEIGHT * attenuation))
						/ INT_ATT_WEIGHT;
	
	return(intScore);	
}

/******************************FindBYIons***********************************
*
*	FindBYIons identifies b and y ions.  
*	The input is as described in the documentation for the function PEFragments, and it 
*	returns a INT_4 containing the
*	value "cleavageSites", which is the number of cleavage sites (amide bonds) that are 
*	defined by a b or y ion.  
*	 The array ionFound is also modified.  
*/
INT_4 FindBYIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				     INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, j, nChargeCount, cChargeCount, bCal, yCal, cleavageSites;
	INT_4 bOrYIon, k, bCalStart, yCalStart;
	INT_4 bIonMass, bIonMassMinErr, bIonMassPlusErr;
	INT_4 yIonMass, yIonMassMinErr, yIonMassPlusErr;
	INT_4 bCount, yCount, bSeries, ySeries, bIon, yIon, massDiff;
	INT_4 precursor, skipOneY, skipOneB;
	INT_4 yCalCorrection = 0;
	INT_4 bCalCorrection = 0;
	char test, twoAAExtension, maxCharge;
	BOOLEAN monoToAvYSwitch = TRUE;	/*Used to recalculate for average masses.*/
	BOOLEAN avToMonoBSwitch = FALSE;	/*Used to recalculate for average masses.*/
	BOOLEAN twoAANTerm;
	REAL_4 currentIonFound;
		
/*	Initialize a few variables.*/
	bCount        = 0;	/*Current number of b ions in a row.*/
	yCount        = 0;	/*Current number of y ions in a row.*/
	bSeries       = 0;	/*The greatest number of b ions in a row for a given sequence.*/
	ySeries       = 0;	/*The greatest number of y ions in a row for a given sequence.*/
	skipOneY      = 1;	/*Allows for one missed y ion in a series w/o resetting to zero.*/
	skipOneB      = 1;	/*Allows for one missed b ion in a series w/o resetting to zero.*/
	cleavageSites = 0;	/*Used to count the number of times a y or b ion of any charge state 
						delineates a sequence. */
	precursor     = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) 
				    / gParam.chargeState;
				    
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}

/*	
	Initialize the b ion starting mass (acetylated, or whatever).
*/
	bCalStart = gParam.modifiedNTerm + 0.5;
	
/*	
	Initialize the y ion starting mass (amidated or unmodified).
*/
	yCalStart = gParam.modifiedCTerm + (2 * gElementMass_x100[HYDROGEN]) + 0.5;

/*	
	Count the number of charged residues.  nChargeCount is one more that the number of 
	charged residues found in an N-terminal fragment. cChargeCount is one more than the number
	of charged residues in a C-terminal fragment.
*/	
	nChargeCount = FindNCharge(sequence, seqLength);
	cChargeCount = 1;

/*	
	Figure out if the N-terminus is a two amino acid extension.
*/
		twoAANTerm = TwoAAExtFinder(sequence, 0);
	
/*
	Here's the big loop, where I step through each position in the sequence.
*/
	for(i = (seqLength - 1); i > 0; i--)	/*Don't do this loop for i = 0 (doesnt make sense).*/
	{
		
/*	
	Initialize some variables for this 'for' loop.
*/
		bOrYIon = 0;	/*If any number of b or y ions of any charge are found, 
						then this equals one.  Otherwise, it stays at zero.*/
		bIon = 0;	/*If a b ion is found this equals one.*/
		yIon = 0;	/*If a y ion is found this equals one.*/
		
/*	
	Figure out if this is a two amino acid extension.
*/
		twoAAExtension = TwoAAExtFinder(sequence, i);
		
/*	
	Calculate the singly charged y ion mass.
*/
	yCal = YCalculator(i, sequence, seqLength, yCalStart, yCalCorrection);
	
/*
	Calculate the singly charged b ion mass.
*/
	bCal = BCalculator(i, sequence, bCalStart, bCalCorrection);
		
/*	
	Readjust the number of charges in the C- and N-terminii.
*/
		if((sequence[i] >= gMonoMass_x100[R] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[R] + gToleranceWide) || 
			(sequence[i] >= gMonoMass_x100[H] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[H] + gToleranceWide) ||
			(sequence[i] >= gMonoMass_x100[K] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[K] + gToleranceWide))
		{
			cChargeCount += 1;
			nChargeCount -= 1;
		}
		else	/*Check to see if its a two amino acid combo that could contain Arg, His, or Lys.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] >= gArgPlus[j] - gToleranceWide 
					&& sequence[i] <= gArgPlus[j] + gToleranceWide) || 
					(sequence[i] >= gHisPlus[j] - gToleranceWide 
					&& sequence[i] <= gHisPlus[j] + gToleranceWide) ||
					(sequence[i] >= gLysPlus[j] - gToleranceWide 
					&& sequence[i] <= gLysPlus[j] + gToleranceWide))
				{
					cChargeCount += 1;
					nChargeCount -= 1;
					break;
				}
			}
		}
		
		for(j = 1; j <= maxCharge; j++) /*Check each charge state up to the parent ion charge.*/
		{

/*	Initialize variables within this j loop.*/
			test = FALSE;	/*Used to test if b, a, or y ions are found before looking 
							for the corresponding losses of ammonia or water.*/
			bIonMass = (bCal + (j * gElementMass_x100[HYDROGEN]) - gElementMass_x100[HYDROGEN]) / j;
				bIonMassMinErr = bIonMass - gToleranceWide;
				bIonMassPlusErr = bIonMass + gToleranceWide;
			yIonMass = (yCal + (j * gElementMass_x100[HYDROGEN]) - gElementMass_x100[HYDROGEN]) / j;
				yIonMassMinErr = yIonMass - gToleranceWide;
				yIonMassPlusErr = yIonMass + gToleranceWide;
			
/*	Search for b ions.*/
			if((bIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				k = fragNum - 1;
				while(fragMOverZ[k] >= bIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= bIonMassPlusErr)
					{
						if(nChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							bOrYIon = 1;	/*A b or y ion has been identified.*/
							test = TRUE;	/*A b ion of charge j has been identified.*/
							bIon = 1;		/*A b ion is present.*/
							massDiff = abs(bIonMass - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
/*	
	If the charge under consideration "j" is equal to the precursor charge, and if that precursor
	charge is greater than one, then the bion intensity is attenuated.  Alternatively, if
	the b ion mass is greater than the precusor, while at the same time the number of basic
	residues in the b ion is less than the number of charges on the precursor, then that is also
	grounds for reducing the influence of the ion.
*/
							if((j == gParam.chargeState && gParam.chargeState > 1) || 
								(bIonMass > precursor && nChargeCount < gParam.chargeState))
							{
								if(gParam.fragmentPattern != 'L')
								{
									ionFound[k] = ionFound[k] * HIGH_MASS_B_ION_MULTIPLIER;
								}
							}
							if(twoAAExtension)
							{
								ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
							}
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}

						}
					}
					k--;
				}
			}

			
/*	Search for the y ion values.*/
			if((yIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				test = FALSE;
				k = fragNum - 1;
				while(fragMOverZ[k] >= yIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= yIonMassPlusErr)
					{
						if(cChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							bOrYIon = 1;	/*A b or y ion has been identified.*/
							test = TRUE;	/*A y ion of charge j has been identified.*/
							yIon = 1;		/*A y ion is present.*/
							massDiff = abs(yIonMass - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							if(j == gParam.chargeState && gParam.chargeState > 1
								&& ((i > 2 && !twoAANTerm) || (i > 1 && twoAANTerm)))
							{
								ionFound[k] = ionFound[k] * HIGH_CHARGE_Y_ION_MULTIPLIER;
							}
							if(twoAAExtension)
							{
								ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
							}
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
						}
					}
					k--;
				}
			}	/*if((yIonMass * j) > ((j-1) * 50000))*/
		}	/*for j*/
		if(bIon)
		{
			bCount++;	/*If there was a b ion, then increment by one.*/
			skipOneB = 1;
		}
		else
		{
			if(skipOneB == 0)
			{
				bCount = 0;	/*Otherwise, reset the counting of b ions to zero.*/
			}
			skipOneB = 0;	/*if the next time there is no b ion, then bCount is reset to zero*/
		}
		if(yIon)
		{
			yCount++;
			skipOneY = 1;
		}
		else
		{
			if(skipOneY == 0)
			{
				yCount = 0;
			}
			skipOneY = 0;	/*if the next time there is no y ion, then bCount is reset to zero*/
		}
		if(bCount > bSeries)
		{
			bSeries = bCount;	/*Don't forget what the longest continuous b series was.*/
		}
		if(yCount > ySeries)
		{
			ySeries = yCount;	/*Don't forget what the longest continuous y series was.*/
		}
		
		if(gParam.fragmentPattern == 'L' && (i < 3 || i > (seqLength - 3)))
		{
			cleavageSites++;	/*For LCQ data, the ends are not well established and should
								not be counted against sequences*/
		}
		else
		{
			cleavageSites += bOrYIon;	/*Count the number of times a b or y ion define a 
										cleavage site.*/
		}
		
	}	/*for i*/
	
	if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q')
	{
		if(ySeries > bSeries)
		{
			cleavageSites = ySeries;
		}
		else
		{
			cleavageSites = bSeries;
		}
	}
	
	return(cleavageSites);
}


/****************************PrintToConsoleAndFile***************************************************
*	This function prints header information plus a list of sequences, scores, ranked according 
*	to the massaged score, which combines the intensity score and x-corr in a single value.  The
*	file that is created is used as the input for the CIDentify database searching program 
*	(a modification of FASTA) written by Alex Taylor.
*/

void	PrintToConsoleAndFile(struct SequenceScore *firstScorePtr, REAL_4 quality, INT_4 length, 
								REAL_4 perfectProbScore)
{
	INT_4 i, j, seqNum, skippedOver = 0, wrongIndex;
	REAL_4 averageWrongXCorrScore, averageWrongIntScore, averageWrongQualityScore, averageWrongProbScore;
	REAL_4 averageWrongComboScore;
	REAL_4 sumOfDiffSquared, diffSquared, stDevXCorrScore, stDevProbScore, stDevIntScore;
	REAL_4 stDevQualityScore, stDevComboScore;
	REAL_4 stDevAboveAvXCorr, stDevAboveAvProb, stDevAboveAvInt;
	REAL_4 stDevAboveAvQuality, stDevAboveAvCombo;
	REAL_4 xCorrRatio, probRatio, intRatio, qualityRatio, comboRatio;
	REAL_4 wrongXCorrScore[WRONG_SEQ_NUM + 1], wrongIntScore[WRONG_SEQ_NUM + 1];
	REAL_4 wrongProbScore[WRONG_SEQ_NUM + 1], wrongQualityScore[WRONG_SEQ_NUM + 1];
	REAL_4 wrongComboScore[WRONG_SEQ_NUM + 1];
	char *peptideString = NULL;
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	INT_4 peptideLength = 0;
	struct SequenceScore *maxPtr;
	FILE *fp;
   	const 	time_t		theTime = (const time_t)time(NULL);
	
        PrintHeaderToFile();

        /* Open the output file for appending.*/
	fp = fopen(gParam.outputFile, "a");
	if(fp == NULL)	/*fopen returns NULL if there's a problem.*/
	{
		printf("Cannot open %s to write the output.\n", gParam.outputFile);
		exit(1);
	}

/*	fprintf(fp, "Spectral Quality = %f\n", quality);
//	fprintf(fp, "Contiguous series of sequence ions defines a sequence of length %2ld\n ", length);
//	fprintf(fp, "High Probscr %f\n ", perfectProbScore);*/
	
	/*fprintf(fp, "Contiguous series of sequence ions defining a sequence of single amino acids %2ld\n ", 
			gSingleAACleavageSites);*/
	fprintf(fp, "\n ");
	fprintf(fp, "\n ");
	fprintf(fp, "\n ");
	
/*        if (gParam.fMonitor) 
//        {
//	  	  	printf("\nMaximum Spectral Quality = %f\n", quality);
//	  	  	printf("\nLongest contiguous series of sequence ions defining a sequence of single amino acids %2ld\n ", length);
//			printf("\nHigh Probscr %f\n ", perfectProbScore);*/
	/*printf("Contiguous series of sequence ions defining a sequence of single amino acids %2ld\n ", 
//			gSingleAACleavageSites);*/
       /* }*/

	
/*	Count the sequences and determine the max seq length.*/
	seqNum = 0;
	maxPtr = firstScorePtr;
	while(maxPtr != NULL)
	{
		seqNum++;
		maxPtr = maxPtr->next;
	}


/*	Set up to print some of the output.*/
	if(gParam.fMonitor)
	{
		printf("\n Sequence                                              Rank  Pr(c)   PevzScr  Quality IntScr X-corr\n");
	}	

	fprintf(fp, "\n Sequence                                              Rank  Pr(c)   PevzScr  Quality IntScr X-corr\n");
	
	if(firstScorePtr == NULL) 
	{
		if(gParam.fMonitor)
		{
			printf("\nNo sequences were found exceeding the specified Pr(c) of %f.\n\n", gParam.outputThreshold);
		}	

		fprintf(fp, "\nNo sequences were found exceeding the specified Pr(c) of %f.\n\n.", gParam.outputThreshold);
		fclose(fp);
		return;
	}
	
/*	Determine if any database sequences looked good, and report their sequences*/

	if (strlen(gParam.databaseSequences) > 0 && gCorrectMass)	/*was a database sequence file opened?*/
	{
		if(gDatabaseSeqCorrect)	/*a database derived sequence looked good*/
		{
			maxPtr = firstScorePtr;
			while(maxPtr != NULL)
			{
				if(maxPtr->databaseSeq == 2)	/*this specific sequence looked good*/
				{
					/*Change peptide[j] to single letter code.*/
					peptideString = NULL;
					peptideLength = 0;
					j = 0;
					while(maxPtr->peptide[j] != 0)
					{
						peptide[j] = maxPtr->peptideSequence[j];
						peptideLength++;
						j++;
					}
					
					peptideString = GetDatabaseSeq(peptide, peptideLength);
					if(peptideString) 
					{
						strcat(peptideString, "  DB");
						if(gParam.fMonitor)
						{
							printf("%-55.55s %2ld   %5.3f   %5.3f    %5.3f   %5.3f  %5.3f\n", 
						  		 peptideString, maxPtr->rank - skippedOver, maxPtr->comboScore, 
								 maxPtr->probScore, maxPtr->quality, maxPtr->intensityScore, 
								 maxPtr->crossDressingScore);
						}
	
						fprintf(fp, "%-55.55s %2ld   %5.3f   %5.3f    %5.3f   %5.3f  %5.3f\n", 
						  		 peptideString, maxPtr->rank - skippedOver, maxPtr->comboScore, 
								 maxPtr->probScore, maxPtr->quality, maxPtr->intensityScore, 
								 maxPtr->crossDressingScore);
								 	 
						free(peptideString);
					}
				}
				maxPtr = maxPtr->next;
			}
		}
	}
	
	for(i = 1; i <= 50; i ++)	/*List the top 50 sequences.*/
	{
		if(i > seqNum)	/*Break out if there are less than 50 sequences in the entire list.*/
		{
			break;
		}
		maxPtr = firstScorePtr;
		while(maxPtr != NULL)
		{	
			if(maxPtr->rank == i && maxPtr->databaseSeq != 2)	
			{
			 	/*Change peptide[j] to single letter code.*/
				peptideString = NULL;
				peptideLength = 0;
				
				j = 0;
				while(maxPtr->peptide[j] != 0)
				{
					peptide[j] = maxPtr->peptideSequence[j];
					peptideLength++;
					j++;
				}
				
				if(maxPtr->databaseSeq)
				{
					peptideString = GetDatabaseSeq(peptide, peptideLength);
				}
				else
				{
					peptideString = PeptideString(peptide, peptideLength);
				}
				
				
								
				
				if(peptideString) 
				{
					if(maxPtr->databaseSeq)
					{
						strcat(peptideString, "  DB");
					}
						 
					if(gParam.fMonitor)
					{
						printf("%-55.55s %2ld   %5.3f   %5.3f    %5.3f   %5.3f  %5.3f\n", 
					  		 peptideString, i - skippedOver, maxPtr->comboScore, 
							 maxPtr->probScore, maxPtr->quality, maxPtr->intensityScore, 
							 maxPtr->crossDressingScore);
					}

					fprintf(fp, "%-55.55s %2ld   %5.3f   %5.3f    %5.3f   %5.3f  %5.3f\n", 
					  		 peptideString, i - skippedOver, maxPtr->comboScore, 
							 maxPtr->probScore, maxPtr->quality, maxPtr->intensityScore, 
							 maxPtr->crossDressingScore);
							 	 
					free(peptideString);
				}
				else
				{
					skippedOver++;
				}
				
				break;
			}
			
			maxPtr = maxPtr->next;
		}
	}
	
	
	
	/*	Print out the sequences w/ x-corr greater than 0.9 that were not already listed above.*/

/*	maxPtr = firstScorePtr;
	while(maxPtr != NULL)
	{
		if(maxPtr->rank <= 50)	//Don't list anything already listed above.
		{
			maxPtr = maxPtr->next;
			continue;
		}
		if(maxPtr->crossDressingScore >= 0.9)
		{
			//Change peptide[j] to single letter code.
			char *peptideString;
			INT_4 peptide[MAX_PEPTIDE_LENGTH];
			INT_4 peptideLength = 0;
			
			j = 0;
			while(maxPtr->peptide[j] != 0)
			{
				peptide[j] = maxPtr->peptideSequence[j];
				peptideLength++;
				j++;
			}
			peptideString = PeptideString(peptide, peptideLength);
			if(maxPtr->databaseSeq)
			{
				strcat(peptideString, "  DB");
			}
			if(peptideString) 
			{
					 
				if(gParam.fMonitor)
				{
					printf("%-55.55s %2ld   %5.3f   %5.3f   %5.3f	  %5.3f\n", 
					  		 peptideString, i - skippedOver, maxPtr->intensityScore, 
							 maxPtr->crossDressingScore, maxPtr->intensityOnlyScore, 
							 maxPtr->quality);
				}

				fprintf(fp, "%-55.55s %2ld   %5.3f   %5.3f   %5.3f	  %5.3f\n",
					  		 peptideString, i - skippedOver, maxPtr->intensityScore, 
							 maxPtr->crossDressingScore, maxPtr->intensityOnlyScore, 
							 maxPtr->quality);
						 	 
				free(peptideString);
			}
		}
		maxPtr = maxPtr->next;
	}
	if(gParam.fragmentPattern == 'Q')
	{
	    fprintf(fp, "\nThe residue 'm' signifies oxidized Met.\n");
            if(gParam.fMonitor)
            {
                printf("\nThe residue 'm' signifies oxidized Met.\n");
            }
	}*/

	/* Print the elapsed search time */
	{
		div_t	theHours;
		div_t	theMin;
		
		theHours = div(gParam.searchTime, 3600);
		theMin   = div(theHours.rem, 60);
		fprintf(fp, "\nSearch time: %2d:", theHours.quot);
		if(gParam.fMonitor)
		{
                    printf("\nSearch time: %2d:", theHours.quot);
		}
		if (theMin.quot < 10)
		{
		    fprintf(fp, "0");
                    if(gParam.fMonitor)
		    {
			printf("0");
		    }
		}
		fprintf(fp, "%d:", theMin.quot);
		if(gParam.fMonitor)
		{
		    printf("%d:", theMin.quot);
		}
                if (theMin.rem < 10) 
		{
			fprintf(fp, "0");
			if(gParam.fMonitor)
		    {
			printf("0");
		    }
		}
		fprintf(fp, "%d\n", theMin.rem);
		if(gParam.fMonitor)
		{
		    printf("%d\n", theMin.rem);
		}
	}
	
/*	
	Here's where I figure out the statistics for the scores for the wrong answers 
	and compare them to the correct mass sequences.
*/
	/*Adjust the number of wrong sequences to avoid counting ones that just had zero scores*/
	wrongIndex = 0;
	for(i = 0; i < gWrongIndex; i++)
	{
		if(	gWrongXCorrScore[i]		!= 0 ||
			gWrongIntScore[i]		!= 0 ||
			gWrongProbScore[i]		!= 0 ||
			gWrongQualityScore[i]	!= 0 ||
			gWrongComboScore[i]		!= 0)
		{
			wrongXCorrScore[wrongIndex]		= gWrongXCorrScore[i];
			wrongIntScore[wrongIndex]		= gWrongIntScore[i];
			wrongProbScore[wrongIndex]		= gWrongProbScore[i];
			wrongQualityScore[wrongIndex]	= gWrongQualityScore[i];
			wrongComboScore[wrongIndex]		= gWrongComboScore[i];
			wrongIndex++;
		}
	}
	/*initialize*/
	for(i = 0; i < WRONG_SEQ_NUM + 1; i++)
	{
		gWrongXCorrScore[i]		= 0;
		gWrongIntScore[i]		= 0;
		gWrongProbScore[i]		= 0;
		gWrongQualityScore[i]	= 0;
		gWrongComboScore[i]		= 0;
	}
	for(i = 0; i < wrongIndex; i++)
	{
		gWrongXCorrScore[i]		= wrongXCorrScore[i];
		gWrongIntScore[i]		= wrongIntScore[i];
		gWrongProbScore[i]		= wrongProbScore[i];
		gWrongQualityScore[i]	= wrongQualityScore[i];
		gWrongComboScore[i]		= wrongComboScore[i];
	}
	gWrongIndex = wrongIndex;
	
	/*If there are enough wrong sequences, then figger out the standard deviations, etc*/
	if(gWrongIndex < 3)
	{
		if(gParam.wrongSeqNum > 0)
		{
			printf("\nNot enough wrong sequences to statistically evaluate\n");
		}
	}
	else
	{
		/*initialize*/
		averageWrongXCorrScore		= 0;
		averageWrongIntScore		= 0;
		averageWrongProbScore		= 0;
		averageWrongQualityScore	= 0;
		averageWrongComboScore		= 0;
		
		/*gWrongIndex-1 (the last in the list) contains the correct sequence score*/
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			averageWrongXCorrScore		+= gWrongXCorrScore[i];
			averageWrongIntScore		+= gWrongIntScore[i];
			averageWrongProbScore		+= gWrongProbScore[i];
			averageWrongQualityScore	+= gWrongQualityScore[i];
			averageWrongComboScore		+= gWrongComboScore[i];
		}
		/*divide by the number of wrong answers to calc the average wrong scores*/
		averageWrongXCorrScore		= averageWrongXCorrScore	/ (gWrongIndex-1);
		averageWrongIntScore		= averageWrongIntScore		/ (gWrongIndex-1);
		averageWrongProbScore		= averageWrongProbScore		/ (gWrongIndex-1);
		averageWrongQualityScore	= averageWrongQualityScore	/ (gWrongIndex-1);
		averageWrongComboScore		= averageWrongComboScore	/ (gWrongIndex-1);
	
/*figure out standard deviations of the wrong cross-correlation scores*/
		sumOfDiffSquared = 0;
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			diffSquared = gWrongXCorrScore[i] - averageWrongXCorrScore;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
		stDevXCorrScore = sumOfDiffSquared / (gWrongIndex - 2);	/*divide by N-1*/
		stDevXCorrScore = sqrt(stDevXCorrScore);
		/*Figure out how many standard deviations the correct mass sequence cross-correlation
		score is above the above the average wrong cross-correlation score*/
		if(stDevXCorrScore != 0)
		{
			stDevAboveAvXCorr = (gWrongXCorrScore[gWrongIndex-1] - averageWrongXCorrScore) 
								/ stDevXCorrScore;
		}
		else
		{
			stDevAboveAvXCorr = 0;	/*avoid a divide by zero*/
		}
		
/*figure out standard deviations of the wrong "probability" scores*/
		sumOfDiffSquared = 0;
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			diffSquared = gWrongProbScore[i] - averageWrongProbScore;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
		stDevProbScore = sumOfDiffSquared / (gWrongIndex - 2);	/*divide by N-1*/
		stDevProbScore = sqrt(stDevProbScore);
		/*Figure out how many standard deviations the correct mass sequence intensity only
		score is above the above the average wrong intensity only score*/
		if(stDevProbScore != 0)
		{
			stDevAboveAvProb = (gWrongProbScore[gWrongIndex-1] - averageWrongProbScore) 
									/ stDevProbScore;
		}
		else
		{
			stDevAboveAvProb = 0;	/*avoid a divide by zero*/
		}
		
/*figure out standard deviations of the wrong intensity scores*/
		sumOfDiffSquared = 0;
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			diffSquared = gWrongIntScore[i] - averageWrongIntScore;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
		stDevIntScore = sumOfDiffSquared / (gWrongIndex - 2);
		stDevIntScore = sqrt(stDevIntScore);
		/*Figure out how many standard deviations the correct mass sequence biased intensity
		score is above the above the average wrong biased intensity score*/
		if(stDevIntScore != 0)
		{
			stDevAboveAvInt = (gWrongIntScore[gWrongIndex-1] - averageWrongIntScore) 
							/ stDevIntScore;
		}
		else
		{	
			stDevAboveAvInt = 0;	/*avoid a divide by zero*/
		}
		
/*figure out standard deviations of the wrong quality scores*/
		sumOfDiffSquared = 0;
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			diffSquared = gWrongQualityScore[i] - averageWrongQualityScore;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
		stDevQualityScore = sumOfDiffSquared / (gWrongIndex - 2);
		stDevQualityScore = sqrt(stDevQualityScore);
		/*Figure out how many standard deviations the correct mass sequence biased intensity
		score is above the above the average wrong biased intensity score*/
		if(stDevQualityScore != 0)
		{
			stDevAboveAvQuality = (gWrongQualityScore[gWrongIndex-1] - averageWrongQualityScore) 
									/ stDevQualityScore;
		}
		else
		{	
			stDevAboveAvQuality = 0;	/*avoid a divide by zero*/
		}
		
/*figure out standard deviations of the wrong combined scores*/
		sumOfDiffSquared = 0;
		for(i = 0; i < gWrongIndex - 1; i++)
		{
			diffSquared = gWrongComboScore[i] - averageWrongComboScore;
			diffSquared = diffSquared * diffSquared;
			sumOfDiffSquared += diffSquared;
		}
		stDevComboScore = sumOfDiffSquared / (gWrongIndex - 2);
		stDevComboScore = sqrt(stDevComboScore);
		/*Figure out how many standard deviations the correct mass sequence biased intensity
		score is above the above the average wrong biased intensity score*/
		if(stDevComboScore != 0)
		{
			stDevAboveAvCombo = (gWrongComboScore[gWrongIndex-1] - averageWrongComboScore) 
									/ stDevComboScore;
		}
		else
		{	
			stDevAboveAvCombo = 0;	/*avoid a divide by zero*/
		}
		
		/*Calculate right/wrong ratios*/
		if(averageWrongXCorrScore != 0)
		{
			xCorrRatio = gWrongXCorrScore[gWrongIndex-1] / averageWrongXCorrScore;
		}
		else
		{
			xCorrRatio = 0;	/*avoid divide by zero*/
		}
		
		if(averageWrongProbScore != 0)
		{
			probRatio = gWrongProbScore[gWrongIndex-1] / averageWrongProbScore;
		}
		else
		{
			probRatio = 0;	/*avoid divide by zero*/
		}
		
		if(averageWrongIntScore != 0)
		{
			intRatio = gWrongIntScore[gWrongIndex-1] / averageWrongIntScore;
		}
		else
		{
			intRatio = 0; /*avoid divide by zero*/
		}
		
		if(averageWrongQualityScore != 0)
		{
			qualityRatio = gWrongQualityScore[gWrongIndex-1] / averageWrongQualityScore;
		}
		else
		{
			qualityRatio = 0; /*avoid divide by zero*/
		}
		
		if(averageWrongComboScore != 0)
		{
			comboRatio = gWrongComboScore[gWrongIndex-1] / averageWrongComboScore;
		}
		else
		{
			comboRatio = 0; /*avoid divide by zero*/
		}
		
		if(gParam.fMonitor)
		{		
			printf("\nStatistics based on %2ld wrong sequences:\n", gWrongIndex-1);
		
            printf("\n                1st ranked  St Deviations   Average Wrong   Correct/Wrong\n");

			printf("xcorr score       %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongXCorrScore[gWrongIndex-1], stDevAboveAvXCorr, averageWrongXCorrScore, xCorrRatio);
			printf("Intensity score   %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongIntScore[gWrongIndex-1], stDevAboveAvInt, averageWrongIntScore, intRatio);
			printf("Quality score     %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongQualityScore[gWrongIndex-1], stDevAboveAvQuality, averageWrongQualityScore, qualityRatio);
			printf("Prob score        %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongProbScore[gWrongIndex-1], stDevAboveAvProb, averageWrongProbScore, probRatio);
			printf("Combined score    %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongComboScore[gWrongIndex-1], stDevAboveAvCombo, averageWrongComboScore, comboRatio);
		}
		
		fprintf(fp, "\nBased on %2ld wrong sequences\n", gWrongIndex);
		
        fprintf(fp, "\n               1st ranked  St Deviations   Average Wrong   Correct/Wrong\n");

		fprintf(fp, "xcorr score      %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongXCorrScore[gWrongIndex-1], stDevAboveAvXCorr, averageWrongXCorrScore, xCorrRatio);
		fprintf(fp, "Intensity score  %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongIntScore[gWrongIndex-1], stDevAboveAvInt, averageWrongIntScore, intRatio);
		fprintf(fp, "Quality score    %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongQualityScore[gWrongIndex-1], stDevAboveAvQuality, averageWrongQualityScore, qualityRatio);
		fprintf(fp, "Prob score       %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongProbScore[gWrongIndex-1], stDevAboveAvProb, averageWrongProbScore, probRatio);
		fprintf(fp, "Combined score   %5.3f          %5.3f          %5.3f          %5.2f\n", 
				gWrongComboScore[gWrongIndex-1], stDevAboveAvCombo, averageWrongComboScore, comboRatio);
		
	}
	fclose(fp);
	return;
}



/*****************************KeepSequence**************************************************************
*
*	If a sequence has an unusually high value for one of the scores, then it is kept regardless of how
*	bad the other scores might be.
*/
BOOLEAN	KeepSequence(struct SequenceScore *currPtr, REAL_4 intscrKeep, REAL_4 xcorrKeep, 
						REAL_4 qualityKeep, REAL_4 probscrKeep)
{
	BOOLEAN	keep = FALSE;
	
	if(currPtr->intensityScore > intscrKeep)
	{
		keep = TRUE;
	}
	if(currPtr->crossDressingScore > xcorrKeep)
	{
		keep = TRUE;
	}
	if(currPtr->quality > qualityKeep)
	{
		keep = TRUE;
	}
	if(currPtr->probScore > probscrKeep)
	{
		keep = TRUE;
	}


	return(keep);
}

/*****************************ComboScore****************************************************************
*
*	ComboScore is derived from an emperical determination of the probability of being wrong given
*	different input scores.  The different probabilities are averaged and given as an output.
*
*/
REAL_4	ComboScore(struct SequenceScore *currPtr)
{
	REAL_4	averageScore, upperLimit, lowerLimit;
	REAL_4	aConstant, bConstant, cConstant;
	REAL_4	probWrong, probRight;
	REAL_4	pevScrWt, qualScrWt, intScrWt, xcorrScrWt, summedWt;

	
	/*Initialize*/
	if(gParam.fragmentPattern == 'Q')
	{
		pevScrWt	= 1;
		qualScrWt	= 1.8;
		intScrWt	= 2;
		xcorrScrWt	= 1;
	}
	else if(gParam.fragmentPattern == 'L')
	{
		pevScrWt	= 1.25;
		qualScrWt	= 0.75;
		intScrWt	= 2;
		xcorrScrWt	= 1.5;
	}
	else
	{
		pevScrWt	= 1;
		qualScrWt	= 1;
		intScrWt	= 1;
		xcorrScrWt	= 1;
	}
	
	summedWt = pevScrWt + qualScrWt + intScrWt + xcorrScrWt;
	if(summedWt == 0)
	{
		printf("Divide by zero in ComboScore");
		exit(1);
	}
	
	/*Calculate the weighted average score*/
	averageScore = currPtr->probScore * pevScrWt + currPtr->intensityScore * intScrWt + 
					currPtr->crossDressingScore * xcorrScrWt + currPtr->quality * qualScrWt;
	averageScore = averageScore / summedWt;
	if(averageScore == 0)
	{
		return(0);	/*return a comboscr of zero*/
	}
	upperLimit			= 0.99;
	lowerLimit			= 0.01;
	
	/*Pr(wrong) = Ax2 + Bx + C*/
	if(gParam.fragmentPattern == 'Q')
	{
		aConstant	= 2.667;	/*Qtof model derived from data set*/
		bConstant	= -4.9901;
		cConstant	= 2.293;
	}
	else if(gParam.fragmentPattern == 'L')
	{
		aConstant	= -1.106;	/*LCQ model derived from data set*/
		bConstant	= -0.607;
		cConstant	= 1.467;
	}
	else
	{
		aConstant	= 0;	/*For other data use an average of the two constants, since no modeling has been done*/
		bConstant	= -2.8;	/*Use a linear model, too*/
		cConstant	= 1.9;
	}
	
	/*Calculate probability of being wrong*/
	probWrong	= aConstant * averageScore * averageScore + bConstant * averageScore + cConstant;
	
	/*Set extreme ends of probability*/
	if(probWrong > upperLimit)
	{
		probWrong = upperLimit;
	}
	else if(probWrong < lowerLimit)
	{
		probWrong = lowerLimit;
	}
	
	
	
	probRight = 1 - probWrong;	/*make it probability of being correct*/
	
	
	
	return(probRight);
}
/*****************************DetermineBestCandidates***************************************************
*
*	Sequences derived from databases are deemed to be correct if either the xcorr, intscr, or probscr
*		values are within 0.95 times the maximum.  If its determined to be correct, then no further scoring
*		occurs.  Otherwise, the scoring proceeds as follows:
*	Sequences with x-corr values less than 0.75 times the max are discarded.
*	Sequences with intensity scores less than 0.75 times the max are discarded.
*	Sequences with quality scores less than max - 0.2 are discarded.
*	Sequences with probscr scores less than max - 0.2 are discarded.
*	A final score is calculated by multiplying quality and probscr, and sequences less than 0.9 of the max
*		are discarded.  
*	If more than five sequences are remaining, then they are ranked by qual x probscr and the top five are
*		kept as the final list of candidates.
*	
*/

struct SequenceScore *DetermineBestCandidates(struct SequenceScore *firstScorePtr)
{
	REAL_4	intscrLimit, intscrBottom, highestIntScore;
	REAL_4	xcorrLimit, xcorrBottom, highestXcorr;
	REAL_4	qualityLimit, qualityBottom, highestQuality;
	REAL_4	probscrLimit, probscrBottom, highestProbscr;
	REAL_4	comboLimit, comboBottom, highestCombo;
	REAL_4	intScore, intOnlyScore, crossDressingScore, stDevErr, calFactor, quality;
	REAL_4	probScore, comboScore;
	REAL_4	intscrKeep, xcorrKeep, qualityKeep, probscrKeep;
	
	INT_4	maxFinalSequences, countTheSeqs, i, cleavageSites;
	INT_4 	sequence[MAX_PEPTIDE_LENGTH], charSequence[MAX_PEPTIDE_LENGTH];
	INT_4	seqLength;

	char	databaseSeq;
	
	BOOLEAN	keep;
	
	struct SequenceScore *currPtr, *massagedSeqListPtr, *previousPtr;
/*
	Initialize variables.
*/
	if(gParam.fragmentPattern == 'Q')
	{
		intscrLimit 		= 0.75;	/*0.75 fraction of highest intscr that is used as threshold*/
		intscrBottom 		= 0.4;	/*0.4 this is the bottom intscr value*/
		intscrKeep			= 0.85; /*keep candidates with intscr greater than this*/
		highestIntScore		= 0;	/*this will hold the max intscr value*/
		
		xcorrLimit 			= 0.75;	/*0.75 fraction of highest xcorr that is used as threshold*/
		xcorrBottom 		= 0.3;	/*0.3 this is the bottom xcorr value*/
		xcorrKeep			= 0.75; /*keep candidates with xcorr greater than this*/
		highestXcorr		= 0;	/*this will hold the max xcorr value*/
		
		qualityLimit 		= 0.35;	/*0.3 this is subtracted from max quality to determine the limit*/
		qualityBottom		= 0.2;	/*0.2 this is the bottom quality value*/
		qualityKeep			= 0.85; /*keep candidates with quality greater than this*/
		highestQuality		= 0;	/*this will hold the max quality value*/
		
		probscrLimit		= 0.2;	/*0.2 this is subtracted from max probscr to determine the limit*/
		probscrBottom		= 0.1;	/*0.1 this is the bottom probscr value*/
		probscrKeep			= 0.75;	/*keep candidates with probscr greater than this*/
		highestProbscr		= 0;	/*this will hold the max probscr value*/
		
		comboLimit 			= 0.8;	/*0.8 fraction of highest probability estimate that is used as threshold*/
		comboBottom			= gParam.outputThreshold;	/*0.2 this is the bottom combo score value*/
		highestCombo		= 0;	/*this will hold the max combo score value*/
	}
	else if(gParam.fragmentPattern == 'L')
	{
		intscrLimit 		= 0.75;	/*0.75 fraction of highest intscr that is used as threshold*/
		intscrBottom 		= 0.45;	/*0.45 this is the bottom intscr value*/
		intscrKeep			= 0.85; /*keep candidates with intscr greater than this*/
		highestIntScore		= 0;	/*this will hold the max intscr value*/
		
		xcorrLimit 			= 0.75;	/*0.75 fraction of highest xcorr that is used as threshold*/
		xcorrBottom 		= 0.35;	/*0.35 this is the bottom xcorr value*/
		xcorrKeep			= 0.85; /*keep candidates with intscr greater than this*/
		highestXcorr		= 0;	/*this will hold the max xcorr value*/
		
		qualityLimit 		= 0.4;	/*0.4 this is subtracted from max quality to determine the limit*/
		qualityBottom		= 0.25;	/*0.25 this is the bottom quality value*/
		qualityKeep			= 0.85; /*keep candidates with intscr greater than this*/
		highestQuality		= 0;	/*this will hold the max quality value*/
		
		probscrLimit		= 0.2;	/*0.2 this is subtracted from max probscr to determine th limit*/
		probscrBottom		= 0.25;	/*0.25 this is the bottom probscr value*/
		probscrKeep			= 0.9;	/*keep candidates with intscr greater than this*/
		highestProbscr		= 0;	/*this will hold the max probscr value*/
		
		comboLimit 			= 0.8;	/*0.8 fraction of highest probability estimate that is used as threshold*/
		comboBottom			= gParam.outputThreshold;	/*0.2 this is the bottom combo score value*/
		highestCombo		= 0;	/*this will hold the max combo score value*/
	}
	else
	{
		intscrLimit 		= 0.75;	/*0.75 fraction of highest intscr that is used as threshold*/
		intscrBottom 		= 0.0;	/*0.0 this is the bottom intscr value*/
		intscrKeep			= 1.0; 	/*keep candidates with intscr greater than this*/
		highestIntScore		= 0;	/*this will hold the max intscr value*/
		
		xcorrLimit 			= 0.75;	/*0.75 fraction of highest xcorr that is used as threshold*/
		xcorrBottom 		= 0.0;	/*0.0 this is the bottom xcorr value*/
		xcorrKeep			= 1.0; 	/*keep candidates with intscr greater than this*/
		highestXcorr		= 0;	/*this will hold the max xcorr value*/
		
		qualityLimit 		= 0.2;	/*0.2 this is subtracted from max quality to determine the limit*/
		qualityBottom		= 0.0;	/*0.0 this is the bottom quality value*/
		qualityKeep			= 1.0; 	/*keep candidates with intscr greater than this*/
		highestQuality		= 0;	/*this will hold the max quality value*/
		
		probscrLimit		= 0.2;	/*0.2 this is subtracted from max probscr to determine th limit*/
		probscrBottom		= 0.0;	/*0.0 this is the bottom probscr value*/
		probscrKeep			= 1.0;	/*keep candidates with intscr greater than this*/
		highestProbscr		= 0;	/*this will hold the max probscr value*/
		
		comboLimit 			= 0.8;	/*0.8 fraction of highest probability estimate that is used as threshold*/
		comboBottom			= gParam.outputThreshold;	/*0.2 this is the bottom combo score value*/
		highestCombo		= 0;	/*this will hold the max combo score value*/
	}
	maxFinalSequences 	= gParam.outputSeqNum;	/*maximum number of sequences in the final list*/
	massagedSeqListPtr 	= NULL;	/*this is the final list of candidate sequences*/
	
	/*Reset bottom limits if a very low output threshold has been selected*/
	if(intscrBottom > gParam.outputThreshold)
	{
		intscrBottom = gParam.outputThreshold;
	}
	if(xcorrBottom > gParam.outputThreshold)
	{
		xcorrBottom = gParam.outputThreshold;
	}
	if(qualityBottom > gParam.outputThreshold)
	{
		qualityBottom = gParam.outputThreshold;
	}
	if(probscrBottom > gParam.outputThreshold)
	{
		probscrBottom = gParam.outputThreshold;
	}
	
	
/*	make sure xcorr, intscr, and quality are less than than one*/
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->crossDressingScore > 1)
		{
			currPtr->crossDressingScore = 1;
		}
		if(currPtr->intensityScore > 1)
		{
			currPtr->intensityScore = 1;
		}
		if(currPtr->quality > 1)
		{
			currPtr->quality = 1;
		}
		currPtr->comboScore = 0;	/*should not have any values in it at this point*/
		currPtr = currPtr->next;
	}
	
/*	Count the sequences and find the max scores.*/
	
	countTheSeqs = 0;
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		countTheSeqs++;
		if(currPtr->intensityScore > highestIntScore)
		{
			highestIntScore = currPtr->intensityScore;
		}
		if(currPtr->crossDressingScore > highestXcorr)
		{
			highestXcorr = currPtr->crossDressingScore;
		}
		if(currPtr->quality > highestQuality)
		{
			highestQuality = currPtr->quality;
		}
		if(currPtr->probScore > highestProbscr)
		{
			highestProbscr = currPtr->probScore;
		}
		currPtr = currPtr->next;
	}
	
/*
*	Determine if any database derived sequences are correct
*/
	
	gDatabaseSeqCorrect = FALSE;
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->databaseSeq)
		{
			if(currPtr->intensityScore >= highestIntScore * 0.95 && 
				currPtr->intensityScore >= intscrBottom)
			{
				gDatabaseSeqCorrect = TRUE;	/*global boolean to report back that there was a database sequence present
											that compared well against the de novo sequences*/
				/*currPtr->comboScore = ComboScore(currPtr);*/
				currPtr->databaseSeq = 2;
			}
			else if(currPtr->crossDressingScore >= highestXcorr * 0.95 && 
				currPtr->crossDressingScore >= xcorrBottom)
			{
				gDatabaseSeqCorrect = TRUE;
				/*currPtr->comboScore = ComboScore(currPtr);*/
				currPtr->databaseSeq = 2;
			}
			else if(currPtr->probScore >= highestProbscr * 0.95 && 
				currPtr->probScore >= probscrBottom)
			{
				gDatabaseSeqCorrect = TRUE;
				/*currPtr->comboScore = ComboScore(currPtr);*/
				currPtr->databaseSeq = 2;
			}
		}
		currPtr = currPtr->next;
	}
	
/*
*	If the database derived sequences are not present, or did not match well, then find the good sequences.
*/

	/*if(!databaseSeqCorrect)
	{*/
		/*Find sequences below xcorr limits and assign zero score values to them.*/
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			/*Does the sequence have unusually high xcorr, intscr, quality or probscr?*/
			keep = KeepSequence(currPtr, intscrKeep, xcorrKeep, qualityKeep, probscrKeep);
			if(!keep)
			{
				if(currPtr->crossDressingScore < highestXcorr * xcorrLimit ||
					currPtr->crossDressingScore < xcorrBottom)
				{
					if((!gDatabaseSeqCorrect) && currPtr->databaseSeq)	/*don't whack the database sequence*/
					{
						currPtr->intensityScore		= 0;
						currPtr->crossDressingScore	= 0;
						currPtr->probScore			= 0;
						currPtr->quality			= 0;
						currPtr->comboScore			= 0;
						currPtr->intensityOnlyScore	= 0;
					}
				}
			}
			currPtr = currPtr->next;
		}
		
		/*Find sequences below intscr limits and assign zero score values to them.*/
		/*But first find the newest high intscr (since it might have been weeded out based on xcorr*/
		highestIntScore = 0;
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->intensityScore > highestIntScore)
			{
				highestIntScore = currPtr->intensityScore;
			}
			currPtr = currPtr->next;
		}
		
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			/*Does the sequence have unusually high xcorr, intscr, quality or probscr?*/
			keep = KeepSequence(currPtr, intscrKeep, xcorrKeep, qualityKeep, probscrKeep);
			if(!keep)
			{
				if(currPtr->intensityScore < highestIntScore * intscrLimit ||
					currPtr->intensityScore < intscrBottom)
				{
					if((!gDatabaseSeqCorrect) && currPtr->databaseSeq)	/*don't whack the database sequence*/
					{
						currPtr->intensityScore		= 0;
						currPtr->crossDressingScore	= 0;
						currPtr->probScore			= 0;
						currPtr->quality			= 0;
						currPtr->comboScore			= 0;
						currPtr->intensityOnlyScore	= 0;
					}
				}
			}
			currPtr = currPtr->next;
		}
		/*Find sequences below quality limits and assign zero score values to them.*/
		/*But first find the newest high quality score*/
		highestQuality = 0;
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->quality > highestQuality)
			{
				highestQuality = currPtr->quality;
			}
			currPtr = currPtr->next;
		}
		
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			/*Does the sequence have unusually high xcorr, intscr, quality or probscr?*/
			keep = KeepSequence(currPtr, intscrKeep, xcorrKeep, qualityKeep, probscrKeep);
			if(!keep)
			{
				if(currPtr->quality < highestQuality - qualityLimit ||
					currPtr->quality < qualityBottom)
				{
					if((!gDatabaseSeqCorrect) && currPtr->databaseSeq)	/*don't whack the database sequence*/
					{
						currPtr->intensityScore		= 0;
						currPtr->crossDressingScore	= 0;
						currPtr->probScore			= 0;
						currPtr->quality			= 0;
						currPtr->comboScore			= 0;
						currPtr->intensityOnlyScore	= 0;
					}
				}
			}
			currPtr = currPtr->next;
		}
		/*Find sequences below probscr limits and assign zero score values to them.*/
		/*But first find the newest high probscr*/
		highestProbscr = 0;
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->probScore > highestProbscr)
			{
				highestProbscr = currPtr->probScore;
			}
			currPtr = currPtr->next;
		}
		
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			/*Does the sequence have unusually high xcorr, intscr, quality or probscr?*/
			keep = KeepSequence(currPtr, intscrKeep, xcorrKeep, qualityKeep, probscrKeep);
			if(!keep)
			{
				if(currPtr->probScore < highestProbscr - probscrLimit ||
					currPtr->probScore < probscrBottom)
				{
					if((!gDatabaseSeqCorrect) && currPtr->databaseSeq)	/*don't whack the database sequence*/
					{
						currPtr->intensityScore		= 0;
						currPtr->crossDressingScore	= 0;
						currPtr->probScore			= 0;
						currPtr->quality			= 0;
						currPtr->comboScore			= 0;
						currPtr->intensityOnlyScore	= 0;
					}
				}
			}
			currPtr = currPtr->next;
		}
		/*Calculate the comboScore*/
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			currPtr->comboScore	= ComboScore(currPtr);	/*returns score of zero if all other scores are zero*/
			currPtr = currPtr->next;
		}
		/*Find highest comboscr value*/
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->comboScore > highestCombo)
			{
				highestCombo = currPtr->comboScore;
			}
			currPtr = currPtr->next;
		}
		/*Find sequences below comboscr limits and assign zero score values to them*/
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->comboScore < highestCombo * comboLimit ||
				currPtr->comboScore < comboBottom)
			{
				if((!gDatabaseSeqCorrect) && currPtr->databaseSeq)	/*don't whack the database sequence*/
				{
					currPtr->intensityScore		= 0;
					currPtr->crossDressingScore	= 0;
					currPtr->probScore			= 0;
					currPtr->quality			= 0;
					currPtr->comboScore			= 0;
					currPtr->intensityOnlyScore	= 0;
				}
			}
			currPtr = currPtr->next;
		}
	/*}*/

/*
*	Store the remaining sequences.
*/
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->comboScore != 0)	/*if there is a comboscore, then save it to the massaged list of seq's;
										if a database-derived sequence was found to be good, then all non-
										database sequences have comboscores of zero and are not saved*/
		{
			intScore			= currPtr->intensityScore;
			intOnlyScore 		= currPtr->intensityOnlyScore;
			crossDressingScore	= currPtr->crossDressingScore;
			stDevErr			= currPtr->stDevErr;
			calFactor 			= currPtr->calFactor;
			quality 			= currPtr->quality;
			probScore 			= currPtr->probScore;
			comboScore 			= currPtr->comboScore;
			cleavageSites 		= currPtr->cleavageSites;
			databaseSeq 		= currPtr->databaseSeq;
			
			/*Load sequence and charSequence arrays.*/
			i = 0;
			while(currPtr->peptide[i] != 0)
			{
				sequence[i] = currPtr->peptide[i];
				charSequence[i] = currPtr->peptideSequence[i];
				i++;
			}
			sequence[i] = 0;
			charSequence[i] = 0;
			seqLength = i;
			
			/*Put into the massaged sequence list.*/
			massagedSeqListPtr = AddToSeqScoreList(massagedSeqListPtr, 
							LoadSeqScoreStruct(intScore, intOnlyScore, sequence, 
												charSequence, seqLength, stDevErr, cleavageSites, 
												calFactor, databaseSeq, crossDressingScore, quality, seqLength,
												probScore, comboScore));
		}
		currPtr = currPtr->next;
	}
	


/*	Rank the sequences.*/
	if(massagedSeqListPtr != NULL)
	{
		SeqComboScoreRanker(massagedSeqListPtr);
	}
	
/*	Eliminate sequences ranked more than maxFinalSequences*/
	if(massagedSeqListPtr != NULL)
	{
		currPtr = massagedSeqListPtr;
		while(currPtr != NULL && currPtr->rank > maxFinalSequences)
		{
			previousPtr = currPtr;
			currPtr = currPtr->next;
			free(previousPtr);
		}
		massagedSeqListPtr = currPtr;	/*found the first sequence ranked less than 6*/
		
		/*now weed out the rest*/
		if(massagedSeqListPtr != NULL && massagedSeqListPtr->next != NULL)
		{
			currPtr = massagedSeqListPtr->next;
			previousPtr = massagedSeqListPtr;
		
			while(currPtr != NULL)
			{
				if(currPtr->rank > maxFinalSequences)
				{
					previousPtr->next = currPtr->next;
					free(currPtr);
					currPtr = previousPtr->next;
				}
				else
				{
					previousPtr = previousPtr->next;
					currPtr = currPtr->next;
				}
			}
		}
	}
	
	/*	Fill in the scores for the top ranked incorrect sequences*/
	if(massagedSeqListPtr == NULL)
	{
		gWrongIndex++;
		return(massagedSeqListPtr);
	}
	currPtr = massagedSeqListPtr;
	while(currPtr != NULL)
	{
		if(currPtr->rank == 1)
		{
			gWrongXCorrScore[gWrongIndex]	= currPtr->crossDressingScore/* * topXCorrScore*/;
			gWrongProbScore[gWrongIndex]	= currPtr->probScore;
			gWrongIntScore[gWrongIndex]		= currPtr->intensityScore;
			gWrongQualityScore[gWrongIndex]	= currPtr->quality;
			gWrongComboScore[gWrongIndex]	= currPtr->comboScore;
			gWrongIndex++;
			break;
		}
		currPtr = currPtr->next;
	}
	
	return(massagedSeqListPtr);
}


/****************************CalcIonFound*****************************************************
*
*	This function uses the global value of gTolerance and an input value that corresponds to the
*	absolute value of the mass difference between a calculated ion m/z and an observed m/z.
*	gTolerance corresponds to a fraction of gParams.fragmentErr (the input fragment ion
*	tolerance).  The function returns a value between zero and one.  If the input mass difference
*	is less than gTolerance, then a one is immediately returned.  If the mass difference is
*	greater than gTolerance, then a value between zero and one is returned, depending on how far
*	off the mass difference is.
*
*
*	New version of CalcIonFound is to use an exponential decay w/ certain boundaries. If the mass
*	difference is less than gToleranceNarrow, then it is assumed to be completely correct.  If it
*	is between gToleranceNarrow and gToleranceWide it is an exponential decay.  Ions in error 
*	greater than gToleranceWide don't even get to this point.
*/

REAL_4	CalcIonFound(REAL_4 currentIonFound, INT_4 massDiff)
{
	REAL_4 ionFound;
/*	REAL_4 range;*/
	
	/*range = gToleranceWide - gToleranceNarrow;
	if(range == 0)
	{
		printf("CalcIonFound: range = 0\n");
		exit(1);
	}*/
	
	if(massDiff <= gToleranceNarrow)
	{
		ionFound = 1;
	}
	else
	{
		/*ionFound = gToleranceWide - massDiff;
		ionFound = ionFound / range;*/
		ionFound = exp(-1 * (massDiff*massDiff) / (2 * gToleranceNarrow * gToleranceNarrow));
	}
	if(ionFound < currentIonFound)
	{
		ionFound = currentIonFound;	/*don't ever replace with a worse value*/
	}
	return(ionFound);
}

/****************************PrintToConsole***************************************************
*	This function prints a list of sequences, scores, ranked according to the intensity score.
*/

void	PrintToConsole(struct SequenceScore *firstScorePtr)
{
	INT_4 i, j, seqNum;
	REAL_4 xcorrNormalizer;
	struct SequenceScore *maxPtr, *currPtr;
	
/*	Find the xcorr normalizer.*/
	xcorrNormalizer = 0;
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->crossDressingScore > xcorrNormalizer)
		{
			xcorrNormalizer = currPtr->crossDressingScore;
		}
		currPtr = currPtr->next;
	}
	xcorrNormalizer = 1;
/*	Set up the screen to print some of the output.*/
/*	printf("\n Rank  X-corr  IntScr  IntOnlyScr Quality  PevScr StDevErr CS  CalFact  Sequence\n");*/
	printf("\n Rank  X-corr  IntScr   Quality  PevScr StDevErr CS  CalFact  Sequence\n");
	
/*	Count the sequences.*/
	seqNum = 0;
	maxPtr = firstScorePtr;
	while(maxPtr != NULL)
	{
		seqNum++;
		maxPtr = maxPtr->next;
	}
		
	for(i = 1; i <= 50 && i <= seqNum; i++)	/*List the top 50 sequences.*/
	{
		maxPtr = firstScorePtr;
		while(maxPtr != NULL)
		{	
			if(maxPtr->rank == i)
			{
				/*Change peptide[j] to single letter code.*/
				char *peptideString;
				INT_4 peptide[MAX_PEPTIDE_LENGTH];
				INT_4 peptideLength = 0;
				
				j = 0;
				while(maxPtr->peptide[j] != 0)
				{
					peptide[j] = maxPtr->peptideSequence[j];
					peptideLength++;
					j++;
				}
				peptideString = PeptideString(peptide, peptideLength);
				if(maxPtr->databaseSeq)
				{
					strcat(peptideString, " ");	/*used to denote this was database sequence*/
				}
				if(peptideString) 
				{
					/*printf(" %3ld   %5.3f   %5.3f   %5.3f      %5.3f    %5.3f  %6.4f  %2ld   %8.6f %s\n", i,
						 maxPtr->crossDressingScore / xcorrNormalizer, maxPtr->intensityScore,
						 maxPtr->intensityOnlyScore, maxPtr->quality, maxPtr->probScore, 
						 maxPtr->stDevErr, maxPtr->cleavageSites, 
						 maxPtr->calFactor, peptideString);*/
					printf(" %3ld   %5.3f   %5.3f    %5.3f    %5.3f  %6.4f  %2ld   %8.6f %s\n", i,
						 maxPtr->crossDressingScore / xcorrNormalizer, maxPtr->intensityScore,
						 maxPtr->quality, maxPtr->probScore, 
						 maxPtr->stDevErr, maxPtr->cleavageSites, 
						 maxPtr->calFactor, peptideString);
					free(peptideString);
				}
			
				break;
			}
			
			maxPtr = maxPtr->next;
		}
	}
		
	/*	Print out the sequences w/ x-corr greater than 0.9 that were not already listed above.*/

	maxPtr = firstScorePtr;
	while(maxPtr != NULL)
	{
		if(maxPtr->rank <= 50)	/*Don't list anything already listed above.*/
		{
			maxPtr = maxPtr->next;
			continue;
		}
		if(maxPtr->crossDressingScore /xcorrNormalizer >= 0.9)
		{
			/*Change peptide[j] to single letter code.*/
			char *peptideString;
			INT_4 peptide[MAX_PEPTIDE_LENGTH];
			INT_4 peptideLength = 0;
			
			j = 0;
			while(maxPtr->peptide[j] != 0)
			{
				peptide[j] = maxPtr->peptideSequence[j];
				peptideLength++;
				j++;
			}
			peptideString = PeptideString(peptide, peptideLength);
			if(maxPtr->databaseSeq)
			{
				strcat(peptideString, "  ");	/*used to denote this was database sequence*/
			}
			if(peptideString) 
			{
					printf(" %3ld   %5.3f   %5.3f   %5.3f      %5.3f    %5.3f  %6.4f  %2ld   %8.6f %s\n", i,
						 maxPtr->crossDressingScore / xcorrNormalizer, maxPtr->intensityScore,
						 maxPtr->intensityOnlyScore, maxPtr->quality, maxPtr->probScore, 
						 maxPtr->stDevErr, maxPtr->cleavageSites, 
						 maxPtr->calFactor, peptideString);
					free(peptideString);
			}
		}
		maxPtr = maxPtr->next;
	}
	
	return;
}

/***********************************AddTagBack**************************************************
*
*	This function is called if a specific sequence tag was used.  It inserts the sequence tag back
*	into the sequence found into the peptide field of firstSequencePtr (linked list of structs of
*	type Sequence).  The value of the field peptideLength is correspondingly adjusted.
*
*/

struct Sequence *AddTagBack(struct Sequence *firstSequencePtr)
{
	INT_4 tagLength, i, j;
	INT_4 tagSequenceMass[MAX_PEPTIDE_LENGTH], peptide[MAX_PEPTIDE_LENGTH];
	struct Sequence *currSeqPtr, *previousPtr, *discardPtr;
	REAL_4 nMass, nMassGroup;
	char test;
	REAL_4 thePeptideMW = gParam.peptideMW;
	
/*	
	First I'll find the length of the sequence tag.
*/		
	tagLength = 0;
	while(gParam.tagSequence[tagLength] != 0)
	{
		tagLength++;
	}

		
/*	
	Next I'll find the N-terminal group mass.
*/
	nMassGroup = gParam.modifiedNTerm;
	
/*	
	I also need to convert the sequence tag from a char array to a list of nominal masses.
*/
	i = 0;
	while(gParam.tagSequence[i] != 0)
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(gGapList[j] != 0)
			{
				if(gSingAA[j] == gParam.tagSequence[i])
				{
					tagSequenceMass[i] = gGapList[j];
				}
			}
		}
		i++;
	}
	
/*	
	Don't forget to readjust the peptideMW.
*/
	for(i = 0; i < tagLength; i++)
	{
		thePeptideMW = thePeptideMW + tagSequenceMass[i];
	}
		
/*	
*	Now I'll start looking at each struct of type Sequence, and modifying the peptide and
*	peptideLength fields.
*/
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		for(i = 0; i < currSeqPtr->peptideLength; i++)	/*Save the old sequence.*/
		{
			peptide[i] = currSeqPtr->peptide[i];
		}
		
		nMass = nMassGroup;	/*Initialize nMass, which serves to keep track of the addition of
							N-terminal amino acids.*/
		
/*	In case the tag sequence is N-terminal, do this first:*/

		test = TRUE;	/*test becomes FALSE if the sequence tag is N-terminal.*/
		if(nMass >= gParam.tagNMass - gToleranceWide)
		{
			if(nMass <= gParam.tagNMass + gToleranceWide)
			{
				test = FALSE;	/*The N-terminal mass is the same as the N-terminal group.
								That is, this sequence tag is an N-terminal sequence.*/
				for(j = 0; j < tagLength; j++)
				{
					currSeqPtr->peptide[j] = tagSequenceMass[j];
				}
				currSeqPtr->peptideLength = (currSeqPtr->peptideLength) + tagLength;
				for(j = tagLength; j < currSeqPtr->peptideLength; j++)
				{
					currSeqPtr->peptide[j] = peptide[j - tagLength];
				}
			}
			else
			{
				currSeqPtr->peptideLength = 0;	/*Set the length to zero as a flag
												that something is wrong w/ the sequence.*/
			}
		}

/*	If the tag sequence is not N-terminal, then the following will occur:*/

		if(test)
		{
			i = 0;	/*Used to follow the peptide length.*/
			while(i < currSeqPtr->peptideLength)
			{
				nMass = nMass + (currSeqPtr->peptide[i]);
				if(nMass >= gParam.tagNMass - gToleranceWide)
				{
					if(nMass <= gParam.tagNMass + gToleranceWide)
					{
						for(j = 0; j < tagLength; j++)
						{
							currSeqPtr->peptide[i + 1 + j] = tagSequenceMass[j];
						}
						for(j = (i + 1); j < currSeqPtr->peptideLength; j++)
						{
							currSeqPtr->peptide[j + tagLength] = peptide[j];
						}
						currSeqPtr->peptideLength = (currSeqPtr->peptideLength) + tagLength;
					}
					else
					{
						currSeqPtr->peptideLength = 0;	/*Set the length to zero as a flag
														that something is wrong w/ the sequence.*/
					}
					break;
				}
				i++;
			}
		}
		currSeqPtr = currSeqPtr->next;
	}
	
	gParam.peptideMW = thePeptideMW;

/*	
	Get rid of any sequences that have a peptideLength of zero.
*/
	while(firstSequencePtr != NULL)	/*Find first sequence that does not have peptidelength of zero*/
	{
		if(firstSequencePtr->peptideLength != 0)
		{
			break;
		}
		discardPtr = firstSequencePtr;
		firstSequencePtr = firstSequencePtr->next;
		free(discardPtr);
	}
	if(firstSequencePtr != NULL)
	{
		previousPtr = firstSequencePtr;
		currSeqPtr = firstSequencePtr->next;
		while(currSeqPtr != NULL)
		{
			if(currSeqPtr->peptideLength == 0)
			{
				previousPtr->next = currSeqPtr->next;
				discardPtr = currSeqPtr;
				currSeqPtr = currSeqPtr->next;
				free(discardPtr);
			}
			else
			{
				currSeqPtr = currSeqPtr->next;
				previousPtr = previousPtr->next;
			}
		}
	}
	

	return(firstSequencePtr);
}

/*****************************SeqIntensityRanker****************************************************
*
*	SeqRanker inputs one pointer to a struct of type SequenceScore - firstScorePtr (points 
*	to the first of the scored sequences).  This function returns nothing and modifies 
*	the "rank" field of this linked list of structs.  The highest ranked sequence score 
*	is ranked "1", etc..  
*/
void SeqIntensityRanker(struct SequenceScore *firstScorePtr)
{
	struct SequenceScore *currPtr, *maxPtr;
	INT_4 i, j, k, seqNum;
	REAL_4 maxScore;
	BOOLEAN test;
	
/*	Count the sequences.*/
	currPtr = firstScorePtr;
	seqNum = 0;
	while(currPtr != NULL)
	{
		seqNum += 1;
		currPtr = currPtr->next;
	}
	
/*	Eliminate any poser database sequences that crept in due to K/Q and the like.*/
	
	currPtr = firstScorePtr;	
	while(currPtr != NULL)
	{
		if(currPtr->databaseSeq)
		{
			for(i = 0; i < gSeqNum; i++)
			{
				test = TRUE;
				j = 0;
				while(gDatabaseSeq[i][j] != 0)
				{
					k = currPtr->peptideSequence[j];
					if(k != I && k != L)
					{
						if(gSingAA[k] !=  gDatabaseSeq[i][j])
						{
							test = FALSE;	/*a mismatch!*/
							break;
						}
					}
					j++;
				}
				if(test)	/*this sequence matches*/
				{
					break;
				}
			}
		
			if(!test)	/*none of the database sequences matched*/
			{
				currPtr->intensityScore = 0;	/*strip away designation as a database sequence 
												and give zero score*/
				currPtr->databaseSeq = 0;
			}
		}
		
		currPtr = currPtr->next;
	}
	
/*	Add 0.05 to the intensityScore for database-derived sequences.  This 0.05 will be taken 
	away down below, so that the score remains unaltered.  However, the ranking will favor
	the database sequences.*/
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->databaseSeq)
		{
			currPtr->intensityScore += 0.05;
		}
		currPtr = currPtr->next;
	}
	
	
/*	Rank the MAX_X_CORR_NUM or seqNum highest scoring sequences (which ever is smaller).*/
	i = 1;
	while(i <= MAX_X_CORR_NUM && i <= seqNum)
	{
		maxPtr = firstScorePtr;
		maxScore = -1000000;
		currPtr = firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->intensityScore > maxScore && currPtr->rank == 0)
			{
				maxPtr = currPtr;
				maxScore = currPtr->intensityScore;
			}
			currPtr = currPtr->next;
		}
		
		maxPtr->rank = i;
		i++;
	}
	
/*	Now that they've been ranked, restore the original intensityScore value.*/
	
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->databaseSeq)
		{
			currPtr->intensityScore -= 0.05;
		}
		currPtr = currPtr->next;
	}
	return;
}

/*****************************SeqComboScoreRanker****************************************************
*
*	SeqRanker inputs one pointer to a struct of type SequenceScore - firstScorePtr (points 
*	to the first of the scored sequences).  This function returns nothing and modifies 
*	the "rank" field of this linked list of structs.  The highest ranked sequence score 
*	is ranked "1", etc..  Ranked according to comboScore field value.
*/
void SeqComboScoreRanker(struct SequenceScore *firstScorePtr)
{
	struct SequenceScore *currPtr, *maxPtr;
	INT_4 i, seqNum;
	REAL_4 maxScore, maxIntScore, maxProbScore;
	
/*	Count the sequences.*/
	currPtr = firstScorePtr;
	seqNum = 0;
	while(currPtr != NULL)
	{
		seqNum += 1;
		currPtr = currPtr->next;
	}
		
/*	Rank the MAX_X_CORR_NUM or seqNum highest scoring sequences (which ever is smaller).*/
	i = 1;
	while(i <= MAX_X_CORR_NUM && i <= seqNum)
	{
		maxPtr 			= firstScorePtr;
		maxScore 		= -1000000;
		maxIntScore 	= -1000000;
		maxProbScore 	= -1000000;
		currPtr 		= firstScorePtr;
		while(currPtr != NULL)
		{
			if(currPtr->comboScore > maxScore && currPtr->rank == 0)	/*rank by comboscore*/
			{
				maxPtr 			= currPtr;
				maxScore 		= currPtr->comboScore;
				maxIntScore		= currPtr->intensityScore;
				maxProbScore 	= currPtr->probScore;
			}
			else if(currPtr->comboScore == maxScore && currPtr->rank == 0)
			{
				if(currPtr->probScore > maxProbScore && currPtr->rank == 0)	/*2nd rank by probScore*/
				{
					maxPtr 			= currPtr;
					maxScore 		= currPtr->comboScore;
					maxIntScore 	= currPtr->intensityScore;
					maxProbScore 	= currPtr->probScore;
				}
				else if(currPtr->probScore == maxProbScore && currPtr->rank == 0)
				{
					if(currPtr->intensityScore > maxIntScore && currPtr->rank == 0)	/*3rd by intScore*/
					{
						maxPtr 			= currPtr;
						maxScore 		= currPtr->comboScore;
						maxIntScore 	= currPtr->intensityScore;
						maxProbScore 	= currPtr->probScore;
					}
				}
			}
			currPtr = currPtr->next;
		}
		
		maxPtr->rank = i;
		i++;
	}
	
	return;
}


/**********************FreeAllSequence************************************
*
* 	Used for freeing memory in a linked list.  
*/
void FreeAllSequence(struct Sequence *currPtr)
{
	struct Sequence *freeMePtr;
	while(currPtr != NULL)
	{
		freeMePtr = currPtr;
		currPtr = currPtr->next;
		free(freeMePtr);
	}
	return;
}

/****************************** FindLowestScore**********************************************
*
*	FindLowestScore is a function that inputs a pointer to a struct of type SequenceScore 
*	that is the first in a linked list of such structs.  It searches through the 
*	intensityScore fields of this list of structs and finds the one with the lowest intensity 
*	score.  It returns a pointer to this struct.
*/
struct SequenceScore *FindLowestScore(struct SequenceScore *currPtr)
{
	struct SequenceScore *lowestPtr;
	
	lowestPtr = currPtr;
	
	while(currPtr != NULL)
	{
		if(currPtr->intensityScore < lowestPtr->intensityScore)
		{
			lowestPtr = currPtr;
		}
		currPtr = currPtr->next;
	}
	
	return(lowestPtr);	
}

/******************************AddToSeqScoreList*********************************
*	AddToSeqScoreList makes a linked list of structs of type SequenceScore.  Input variables 
*	are two pointers
*	to structs of type SequenceScore.  The first pointer is firstScorePtr, which is initially 
*	set to NULL, but
*	is given the address of the first struct in the linked list and remains unchanged 
*	thereafter.  The second
*	pointer points to the new struct to be added.  This second pointer already has the 
*	'next' field set to NULL.
*/
struct SequenceScore *AddToSeqScoreList(struct SequenceScore *firstPtr, 
											struct SequenceScore *currPtr)
{
	static struct SequenceScore *lastPtr;
	
	if(firstPtr == NULL)
	{
		firstPtr = currPtr;
	}
	else
	{
		lastPtr->next = currPtr;
	}
	lastPtr = currPtr;
	
	return(firstPtr);
}

/******************************LoadSeqScoreStruct********************************
*
*	LoadSeqScoreStruct inputs "intScore" (the intensity-based score for the sequence), 
*	"sequence" (the nominal mass of the extensions), and "seqLength" (the sequence length).  
*	It returns a pointer to a struct of type SequenceScore, which has had its peptide 
*	and intensityScore fields loaded using the input values above.
*	The remaining fields (crossDressingScore, rank, and next) are NULLed.
*/
struct SequenceScore *LoadSeqScoreStruct(REAL_4 intScore, REAL_4 intOnlyScore,
										INT_4 *sequence, INT_4 *charSequence, INT_4 seqLength,
										REAL_4 stDevErr, INT_4 cleavageSites, REAL_4 calFactor,
										char databaseSeq, REAL_4 normalXCorScore, REAL_4 quality,
										REAL_4 length, REAL_4 probScore, REAL_4 comboScore)
{
	struct SequenceScore *currPtr;
	INT_4 i;
	
	currPtr = (struct SequenceScore *) malloc(sizeof(struct SequenceScore));
	if(currPtr == NULL)
	{
		printf("SequenceScore: Out of memory");
		exit(1);
	}
	
	for(i = 0; i < seqLength; i++)
	{
		currPtr->peptide[i] = sequence[i];
		currPtr->peptideSequence[i] = charSequence[i];
	}
	currPtr->peptide[seqLength] = 0;
	
	currPtr->peptideSequence[seqLength] = 0;
	
	currPtr->intensityScore = intScore;
	
	currPtr->calFactor = calFactor;
	
	currPtr->cleavageSites = cleavageSites;
	
	currPtr->stDevErr = stDevErr;
	
	currPtr->databaseSeq = databaseSeq;
	
	currPtr->probScore = probScore;
	
	currPtr->intensityOnlyScore = intOnlyScore;
	
	currPtr->crossDressingScore = normalXCorScore;
	
	currPtr->comboScore = comboScore;
	
	currPtr->rank = 0;
	
	currPtr->length = length;
	
	currPtr->quality = quality;
	
	currPtr->next = NULL;
	
	return(currPtr);	
}


/******************************IntensityOnlyScorer*********************************
*
*	IntensityScorer inputs fragIntensity (the ion intensities), ionFound (the array that 
*	is indexed the same as fragIntensity and contains 1 for ions that have been identified 
*	and 0 for those that were not), and fragNum (the number of fragment
*	ions in the CID data).  This function returns a REAL_4 value (ranging from zero to one) 
*	corresponding to the fraction of the ion current that can be identified times as one
*	of the known ion types.
*/
REAL_4 IntensityOnlyScorer(INT_4 *fragIntensity, REAL_4 *ionFound, 
						INT_4 fragNum, INT_4 intensityTotal)
{
	REAL_4 intOnlyScore = 0;
	INT_4 i;
	
	
/*	Add up the intensity that has been identified, and count the ions.*/
	
	for(i = 0; i < fragNum; i++)
	{
		if(ionFound[i] != 0)
		{
			intOnlyScore += (REAL_4)fragIntensity[i];
		}
	}	
	if(intensityTotal == 0)
	{
		printf("IntensityOnlyScorer:  intensityTotal = 0\n");
		exit(1);
	}
	intOnlyScore = intOnlyScore / intensityTotal;
	return(intOnlyScore);	
}


/******************************IntensityScorer*********************************
*
*	IntensityScorer inputs fragIntensity (the ion intensities), ionFound (the array that 
*	is indexed the same as fragIntensity and contains 1 for ions that have been identified 
*	and 0 for those that were not), cleavageSites (the number of times either a b or y ion 
*	of any charge was found to delineate the sequence), and fragNum (the number of fragment
*	ions in the CID data).  This function returns a REAL_4 value (ranging from zero to one) 
*	corresponding to the fraction of the ion current that can be identified times a multiplier 
*	that reflects the idea that correct sequences are usually delineated by series of either 
*	b or y ions.
*/
REAL_4 IntensityScorer(INT_4 *fragIntensity, REAL_4 *ionFound, INT_4 cleavageSites, 
						INT_4 fragNum, INT_4 seqLength, INT_4 intensityTotal)
{
	REAL_4 intScore = 0;
	REAL_4 numScore = 0;
	REAL_4 attenuation, avPeaksPerResidue, avResidueNum, peaksPerResidue, realPerAvRatio;
	INT_4 i;
	INT_4 intNum = 0;
	

/*	
*	Initialize attenuation, which is a fractional multiplier that reflects the number of 
*	times either b or y ions delineate
*	the proposed sequence.
*/
	if(seqLength - 1 == 0)
	{
		printf("IntensityScorer: seqLength - 1 = 0\n");
		exit(1);
	}
	attenuation = (REAL_4)cleavageSites / ((REAL_4)(seqLength - 1));
	if(attenuation > 1)
	{
		attenuation = 1;	/*make sure this never exceeds one*/
	}
	if(attenuation < 0)
	{
		attenuation = 0;
	}
	/*attenuation = attenuation / (seqLength - 1);*/
	
/*	Add up the intensity that has been identified, and count the ions.*/
	
	for(i = 0; i < fragNum; i++)
	{
		if(ionFound[i] != 0)
		{
			intScore += (REAL_4)fragIntensity[i] * ionFound[i];
			numScore += ionFound[i];
			/*numScore++;*/
			intNum++;
		}
	}
	
/*	Figure out the number of peaks per average number of residues.*/
	avResidueNum = gParam.peptideMW / (REAL_4)gAvResidueMass;
	if(avResidueNum == 0)
	{
		printf("IntensityScorer: avResidueNum = 0\n");
		exit(1);
	}
	avPeaksPerResidue = (REAL_4)intNum / (REAL_4)avResidueNum;

/*	Figure out the number of peaks per actual residue.*/
	if(seqLength == 0)
	{
		printf("IntensityScorer: seqLength = 0\n");
		exit(1);
	}
	peaksPerResidue = (REAL_4)intNum / (REAL_4)seqLength;
	
/*	Now determine the ratio of peaksPerResidue to avPeaksPerResidue.*/
	if(avPeaksPerResidue == 0)
	{
		return(0);	/*avPeaksPerResidue is 0 if no ions match, so return score of 0*/
	}
	realPerAvRatio = (REAL_4)peaksPerResidue / (REAL_4)avPeaksPerResidue;
	if(realPerAvRatio > 1)
	{
		realPerAvRatio = 1;
	}
	if(realPerAvRatio < 0)
	{
		realPerAvRatio = 0;
	}

/*	Percent of ion current accounted for.*/
	if(intensityTotal == 0)	/*Avoid divide by zero.*/
	{
		printf("IntensityScorer:  intensityTotal = 0\n");
		exit(1);
	}
	intScore = (REAL_4)intScore / (REAL_4)intensityTotal;
	if(intScore > 1)
	{
		intScore = 1;
	}
	if(intScore < 0)
	{
		intScore = 0;
	}
	
/*	Percent of the number of ions accounted for regardless of their intensity.*/
	if(fragNum == 0)	/*Avoid divide by zero.*/
	{
		printf("IntensityScorer:  fragNum = 0\n");
		exit(1);
	}
	numScore = (REAL_4)numScore / (REAL_4)fragNum;
	if(numScore > 1)
	{
		numScore = 1;
	}
	if(numScore < 0)
	{
		numScore = 0;
	}
	
/*	Put the numScore, intScore, attenuation, and realPerAvRatio all together.*/
	intScore = ((INTENSITY_WEIGHT * intScore) + (PEAKS_WEIGHT * realPerAvRatio)
				+ (ATTENUATION_WEIGHT * attenuation) + (NUMBER_WEIGHT * numScore))
				/ ATT_INT_PEAKS_NUM;
	/*intScore = intScore * intScore * attenuation * numScore;*/
	
	return(intScore);	
}



/************************* ScoreLowMassIons ***********************************************
*
*	ScoreLowMassIons inputs the same information as the function PEFragments (plus the array 
*	"lowMassIons"), and  returns nothing.  Only the array "ionFound" is modified.  
*	This function determines which amino acids are present in the sequence, and then uses
*	the array lowMassIons for identification of the low mass amino acid - specific ions.  
*	Only singly charged ions are considered.
*/
void ScoreLowMassIons(REAL_4 *ionFound, INT_4 *fragMOverZ, INT_4 *sequence, 
					  INT_4 seqLength, INT_4 lowMassIons[][3], INT_4 *ionType)
{
	INT_4 i, j, k, m, n, p;
	INT_4 lowMassPlusErr, lowMassMinErr, massDiff;
	REAL_4 currentIonFound;
	char test;
	
	

/*	Using the array "sequence", identify the low mass ions.*/
	for(i = 0; i < seqLength; i++)
	{
		for(j = 0; j < 3; j++)
		{
			for(m = 0; m < gAminoAcidNumber; m++)
			{
				if((sequence[i] <= gMonoMass_x100[m] + gToleranceWide) && 
					(sequence[i] >= gMonoMass_x100[m] - gToleranceWide))
				{
					if(lowMassIons[m][j] != 0)
					{
						k = 0;
						lowMassPlusErr = lowMassIons[m][j] + gToleranceWide;
						lowMassMinErr = lowMassIons[m][j] - gToleranceWide;
						while(fragMOverZ[k] <= lowMassPlusErr)
						{
							if(fragMOverZ[k] >= lowMassMinErr)
							{
								currentIonFound = ionFound[k];
								massDiff = abs(lowMassIons[m][j] - fragMOverZ[k]);
								ionFound[k] = CalcIonFound(ionFound[k], massDiff);
								if(currentIonFound > ionFound[k])
								{
									ionFound[k] = currentIonFound;
								}
								else
								{
									ionType[k] = 14;	/*low mass ions*/
								}
							}
							k++;
						}
					}
				}
			}
		}
	}
	
/*	Check for two aa extensions */
	for(i = 0; i < seqLength; i++)
	{
		test = TRUE;	/*yes it is a two aa extension*/
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if((sequence[i] <= gMonoMass_x100[j] + gToleranceWide) && 
					(sequence[i] >= gMonoMass_x100[j] - gToleranceWide))
			{
				test = FALSE;	/*no, its not a two aa extension*/
			}
		}
		if(test)
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gGapList[j] != 0)
				{
					for(m = 0; m <  gAminoAcidNumber; m++)
					{
						if(gGapList[m] != 0)
						{
							if((sequence[i] <= gGapList[j] + gGapList[m] + gToleranceWide) &&
								(sequence[i] >= gGapList[j] + gGapList[m] - gToleranceWide))
							{
								for(n = 0; n < 3; n++)
								{
									if(lowMassIons[j][n] != 0)
									{
										k = 0;
										lowMassPlusErr = lowMassIons[j][n] + gToleranceWide;
										lowMassMinErr = lowMassIons[j][n] - gToleranceWide;
										while(fragMOverZ[k] <= lowMassPlusErr)
										{
											if(fragMOverZ[k] >= lowMassMinErr)
											{
												currentIonFound = ionFound[k];
												massDiff = abs(lowMassIons[j][n] - fragMOverZ[k]);
												ionFound[k] = CalcIonFound(ionFound[k], massDiff);
												if(currentIonFound > ionFound[k])
												{
													ionFound[k] = currentIonFound;
												}
												else
												{
													ionType[k] = 14;	/*low mass ions*/
												}
											}
											k++;
										}
									}
								}
								for(n = 0; n < 3; n++)
								{
									if(lowMassIons[m][n] != 0)
									{
										k = 0;
										lowMassPlusErr = lowMassIons[m][n] + gToleranceWide;
										lowMassMinErr = lowMassIons[m][n] - gToleranceWide;
										while(fragMOverZ[k] <= lowMassPlusErr)
										{
											if(fragMOverZ[k] >= lowMassMinErr)
											{
												currentIonFound = ionFound[k];
												massDiff = abs(lowMassIons[m][n] - fragMOverZ[k]);
												ionFound[k] = CalcIonFound(ionFound[k], massDiff);
												if(currentIonFound > ionFound[k])
												{
													ionFound[k] = currentIonFound;
												}
												else
												{
													ionType[k] = 14;	/*low mass ions*/
												}
											}
											k++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
/*	Check for three aa extensions */
	for(i = 0; i < seqLength; i++)
	{
		test = TRUE;	/*yes it is a two or three aa extension*/
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if((sequence[i] <= gMonoMass_x100[j] + gToleranceWide) && 
					(sequence[i] >= gMonoMass_x100[j] - gToleranceWide))
			{
				test = FALSE;	/*no, its not a two aa extension*/
			}
		}
		if(test)
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gGapList[j] != 0)
				{
					for(m = 0; m <  gAminoAcidNumber; m++)
					{
						if(gGapList[m] != 0)
						{
							for(p = 0; p < gAminoAcidNumber; p++)
							{
								if(gGapList[p] != 0)
								{
									if((sequence[i] <= gGapList[j] + gGapList[m] + gGapList[p] + gToleranceWide) &&
										(sequence[i] >= gGapList[j] + gGapList[m] + gGapList[p] - gToleranceWide))
									{
										for(n = 0; n < 3; n++)
										{
											if(lowMassIons[j][n] != 0)
											{
												k = 0;
												lowMassPlusErr = lowMassIons[j][n] + gToleranceWide;
												lowMassMinErr = lowMassIons[j][n] - gToleranceWide;
												while(fragMOverZ[k] <= lowMassPlusErr)
												{
													if(fragMOverZ[k] >= lowMassMinErr)
													{
														currentIonFound = ionFound[k];
														massDiff = abs(lowMassIons[j][n] - fragMOverZ[k]);
														ionFound[k] = CalcIonFound(ionFound[k], massDiff);
														if(currentIonFound > ionFound[k])
														{
															ionFound[k] = currentIonFound;
														}
														else
														{
															ionType[k] = 14;	/*low mass ions*/
														}
													}
													k++;
												}
											}
										}
										for(n = 0; n < 3; n++)
										{
											if(lowMassIons[m][n] != 0)
											{
												k = 0;
												lowMassPlusErr = lowMassIons[m][n] + gToleranceWide;
												lowMassMinErr = lowMassIons[m][n] - gToleranceWide;
												while(fragMOverZ[k] <= lowMassPlusErr)
												{
													if(fragMOverZ[k] >= lowMassMinErr)
													{
														currentIonFound = ionFound[k];
														massDiff = abs(lowMassIons[m][n] - fragMOverZ[k]);
														ionFound[k] = CalcIonFound(ionFound[k], massDiff);
														if(currentIonFound > ionFound[k])
														{
															ionFound[k] = currentIonFound;
														}
														else
														{
															ionType[k] = 14;	/*low mass ions*/
														}
													}
													k++;
												}
											}
										}
										for(n = 0; n < 3; n++)
										{
											if(lowMassIons[p][n] != 0)
											{
												k = 0;
												lowMassPlusErr = lowMassIons[p][n] + gToleranceWide;
												lowMassMinErr = lowMassIons[p][n] - gToleranceWide;
												while(fragMOverZ[k] <= lowMassPlusErr)
												{
													if(fragMOverZ[k] >= lowMassMinErr)
													{
														currentIonFound = ionFound[k];
														massDiff = abs(lowMassIons[p][n] - fragMOverZ[k]);
														ionFound[k] = CalcIonFound(ionFound[k], massDiff);
														if(currentIonFound > ionFound[k])
														{
															ionFound[k] = currentIonFound;
														}
														else
														{
															ionType[k] = 14;	/*low mass ions*/
														}
													}
													k++;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
		
	return;
}

/****************************** ProInThirdPosition *******************************
*
*	If there is a proline in the third position for LCQ data, then a value of -1
*	is assigned as the score, which is used later to prevent tossing this sequence
*	out.
*/

REAL_4 ProInThirdPosition(REAL_4 oldScore, INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, k, aaCount;
	REAL_4 score;
	char test;
	
	aaCount = 0;
	for(k = 0; k < seqLength; k++)
	{
		if(k < 4)
		{
			test = TRUE;
			for(i = 0; i < gAminoAcidNumber; i++)
			{
				if(sequence[k] <= gGapList[i] + gToleranceWide &&
					sequence[k] >= gGapList[i] - gToleranceWide)
				{
					test = FALSE;
				}
			}
			if(test)
			{
				aaCount++;
				aaCount++;
			}
			else
			{
				aaCount++;
			}
			if(aaCount >= 3)
				break;
		}
	}
	
	if(aaCount == 3)
	{
		if(sequence[k] <= gGapList[P] + gToleranceWide &&
			sequence[k] >= gGapList[P] - gToleranceWide)
		{
			score = -1;
		}
		else
		{
			score = oldScore;
		}
	}
	else
	{
		score = oldScore;
	}
	
	return(score);
}

/****************************** AssignHighMZScore **********************************
*
*
*
*/
REAL_4 AssignHighMZScore(INT_4 highMZNum, INT_4 *highMZFrags, INT_4 *highMZInts, 
						REAL_4 totalIntensity, INT_4 *sequence, INT_4 seqLength)
{
	INT_4 bCalStart, yCalStart, nChargeCount, cChargeCount;
	INT_4 i, j, k, m, bIonMass, bIonMassMinErr, bIonMassPlusErr;
	INT_4 bMinW, bMinWMinErr, bMinA, bMinAPlusErr;
	INT_4 yMinW, yMinWMinErr, yMinA, yMinAPlusErr;
	INT_4 yIonMass, yIonMassMinErr, yIonMassPlusErr;
	INT_4 yCal, bCal;
	INT_4 yCalCorrection = 0;
	INT_4 bCalCorrection = 0;
	REAL_4 score = 0;
	char test, test1, maxCharge;
	REAL_4 *ionFound;
	
	ionFound = (float *) malloc(MAX_ION_NUM * sizeof(REAL_4));
	if(ionFound == NULL)
	{
		printf("AssignHighMZScore:  Out of memory.");
		exit(1);
	}
	/*Initialize ionFound*/
	for(j = 0; j < MAX_ION_NUM; j++)
	{
		ionFound[j] = 0;
	}
	
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}

/*	
	Initialize the starting b ion mass (acetylated, etc).
*/
	bCalStart = gParam.modifiedNTerm + 0.5;
	
/*	
	Initialize the starting mass for y ions.
*/
	yCalStart = gParam.modifiedCTerm + (2 * gElementMass_x100[HYDROGEN]) + 0.5;	
	
/*	
	Count the number of charged residues.  nChargeCount is one more than the number of 
	charged residues found in an N-terminal fragment.  cChargeCount is one more than the 
	number of charged residues found in a C-terminal fragment.
*/
	nChargeCount = FindNCharge(sequence, seqLength);
	cChargeCount = 1;
	
	/*Step through the sequence.*/
	for(i = (seqLength - 1); i > 0; i--)	/*Don't do this loop for i = 0 (doesnt make sense).*/
	{

/*
	Calculate the singly charged y ion mass.
*/
		yCal = YCalculator(i, sequence, seqLength, yCalStart, yCalCorrection);	
			
/*
	Calculate the singly charged b ion mass.
*/
		bCal = BCalculator(i, sequence, bCalStart, bCalCorrection);		
		
/*	
	Readjust the number of charges in the C- and N-terminii.
*/
		if((sequence[i] >= gMonoMass_x100[R] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[R] + gToleranceWide) || 
			(sequence[i] >= gMonoMass_x100[H] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[H] + gToleranceWide) ||
			(sequence[i] >= gMonoMass_x100[K] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[K] + gToleranceWide))
		{
			cChargeCount += 1;
			nChargeCount -= 1;
		}
		else	/*Check to see if its a two amino acid combo that could contain Arg, His, or Lys.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] >= gArgPlus[j] - gToleranceWide 
					&& sequence[i] <= gArgPlus[j] + gToleranceWide) || 
					(sequence[i] >= gHisPlus[j] - gToleranceWide 
					&& sequence[i] <= gHisPlus[j] + gToleranceWide) ||
					(sequence[i] >= gLysPlus[j] - gToleranceWide 
					&& sequence[i] <= gLysPlus[j] + gToleranceWide))
				{
					cChargeCount += 1;
					nChargeCount -= 1;
					break;
				}
			}
		}
		
/*
	Check each charge state up to the parent ion charge.
*/
		for(j = 1; j <= maxCharge; j++) 
		{

/*	
	Initialize variables within this j loop.
*/
			test = FALSE;	/*Used to test if b, a, or y ions are found before looking 
							for the corresponding losses of ammonia or water.*/
			bIonMass = (bCal + (j * gElementMass_x100[HYDROGEN]) - 
						gElementMass_x100[HYDROGEN]) / j;
				bIonMassMinErr = bIonMass - gToleranceWide;
				bIonMassPlusErr = bIonMass + gToleranceWide;
			yIonMass = (yCal + (j * gElementMass_x100[HYDROGEN]) 
						- gElementMass_x100[HYDROGEN]) / j;
				yIonMassMinErr = yIonMass - gToleranceWide;
				yIonMassPlusErr = yIonMass + gToleranceWide;
			
/*	Search for b ions.*/
			if((bIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				k = highMZNum - 1;
				while(highMZFrags[k] >= bIonMassMinErr && k >= 0)
				{
					if(highMZFrags[k] <= bIonMassPlusErr)
					{
						if(nChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							test = TRUE;	/*A b ion of charge j has been identified.*/
							ionFound[k] = 1;
						}
					}
					k--;
				}

/*	
*	Search for b minus ammonia or water.  The index value of k is carried over from the 
*	b ion search, since
*	the while loop breaks out once the mass value is less than the calculated b value 
*	(ie, k is close to 
*	b minus ammonia or water.
*/
				test1 = FALSE;
				if(!test)
				{
					if(sequence[0] <= gMonoMass_x100[K] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[K] - gParam.fragmentErr
						&& gParam.fragmentErr > 0.04 * gMultiplier)
					{
						test1 = TRUE;
					}
					if(sequence[0] <= gMonoMass_x100[Q] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[Q] - gParam.fragmentErr)
					{
						test1 = TRUE;
					}
					if(!test1)
					{
						for(m = 0; m < gAminoAcidNumber; m++)
						{
							if(sequence[0] <= gGlnPlus[m] + gParam.fragmentErr &&
								sequence[0] >= gGlnPlus[m] - gParam.fragmentErr)
							{
								test1 = TRUE;
								break;
							}
							if(sequence[0] + sequence[1] <= gGlnPlus[m] + gParam.fragmentErr &&
								sequence[0] + sequence[1] >= gGlnPlus[m] - gParam.fragmentErr)
							{
								test1 = TRUE;
								break;
							}
						}
					}
				}
				
				if(!test1)
				{
					if(sequence[0] <= gMonoMass_x100[E] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[E] - gParam.fragmentErr)
					{
						test1 = TRUE;
					}
					
					if(!test1)
					{
						for(m = 0; m < gAminoAcidNumber; m++)
						{
							if(sequence[0] <= gGluPlus[m] + gParam.fragmentErr &&
								sequence[0] >= gGluPlus[m] - gParam.fragmentErr)
							{
								test1 = TRUE;
								break;
							}
							if(sequence[0] + sequence[1] <= gGluPlus[m] + gParam.fragmentErr &&
								sequence[0] + sequence[1] >= gGluPlus[m] - gParam.fragmentErr)
							{
								test1 = TRUE;
								break;
							}
						}
					}
				}
				
				if(test || test1)
				{
				
/*	Calculate the b minus water and b minus ammonia values.*/
					bMinW = (bCal - gWater + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						bMinWMinErr = bMinW - gToleranceWide;
					bMinA = (bCal - gAmmonia + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						bMinAPlusErr = bMinA + gToleranceWide;
					
					if(nChargeCount >= j)
					{
						while(highMZFrags[k] >= bMinWMinErr && k >= 0)
						{
							if(highMZFrags[k] <= bMinAPlusErr)
							{
								if(test)
								{
									ionFound[k] = 1;	/*full count if b ion is present*/
								}
								else
								{
									ionFound[k] = 0.75;	/*partial count if Nterm QorE*/
								}
							}
							k--;
						}
					}
				}
			}	/*if((bIonMass * j) > ((j-1) * 50000))*/
			
/*	Search for the y ion values.*/
			if((yIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				test = FALSE;
				k = highMZNum - 1;
				while(highMZFrags[k] >= yIonMassMinErr && k >= 0)
				{
					if(highMZFrags[k] <= yIonMassPlusErr)
					{
						if(cChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							test = TRUE;	/*A y ion of charge j has been identified.*/
							ionFound[k] = 1;
						}
					}
					k--;
				}

/*	
*	Search for y minus ammonia or water.  The index value of k is carried over from the 
*	y ion search, since
*	the while loop breaks out once the mass value is less than the calculated y value 
*	(ie, k is close to 
*	y minus ammonia or water.
*/
				if(test)
				{
				
/*	Calculate the y minus water and y minus ammonia values.*/
					yMinW = (yCal - gWater + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						yMinWMinErr = yMinW - gToleranceWide;
					yMinA = (yCal - gAmmonia + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						yMinAPlusErr = yMinA + gToleranceWide;
					if(cChargeCount >= j)
					{
						while(highMZFrags[k] >= yMinWMinErr && k >= 0)
						{
							if(highMZFrags[k] <= yMinAPlusErr)
							{
								ionFound[k] = 1;
							}
							k--;
						}
					}
				}
			}	/*if((yIonMass * j) > ((j-1) * 50000))*/
		}	/*for j*/
	}	/*for i*/
	
	score = 0;
	for(i = 0; i < highMZNum; i++)
	{
		if(ionFound[i] != 0)
		{
			score += highMZInts[i] * ionFound[i];
		}
	}
	
	if(totalIntensity == 0)
	{
		printf("AssignHighMZScore: totalIntensity = 0\n");
		exit(1);
	}
	score = score / totalIntensity;

	free(ionFound);
	
	return(score);
}

/******************************FindABYIons***********************************
*
*	FindABYIons identifies a, b, and y ions, plus losses of ammonia and water from these 
*	three ion types.  
*	The input is as described in the documentation for the function PEFragments, and it 
*	returns a INT_4 containing the
*	value "cleavageSites", which is the number of cleavage sites (amide bonds) that are 
*	defined by a b or y ion.  
*	 The array ionFound is also modified.  
*/
INT_4 FindABYIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				     INT_4 *sequence, INT_4 seqLength, char argPresent,
				     REAL_4 *yFound, REAL_4 *bFound, REAL_8 *byError, INT_4 *ionType)

{
	INT_4 i, j, nChargeCount, cChargeCount, bCal, yCal, cleavageSites;
	INT_4 bOrYIon, k, aIonMass, aIonMassMinErr, aIonMassPlusErr;
	INT_4 aMinW, aMinWMinErr, aMinA, aMinAPlusErr, bCalStart, yCalStart;
	INT_4 bIonMass, bIonMassMinErr, bIonMassPlusErr;
	INT_4 bMinW, bMinA, bMinWMinErr, bMinAPlusErr;
	INT_4 yIonMass, yIonMassMinErr, yIonMassPlusErr;
	INT_4 yMinW, yMinWMinErr, yMinA, yMinAPlusErr;
	INT_4 bMin64, lossOf64, bMin64MinErr, bMin64PlusErr;
	INT_4 yMin64, yMin64MinErr, yMin64PlusErr;
	INT_4 bCount, yCount, bSeries, ySeries, bIon, yIon, massDiff, massDiffW, massDiffA;
	INT_4 bCountStringent, yCountStringent, bSeriesStringent, ySeriesStringent;
	INT_4 precursor, m, skipOneY, skipOneB;
	INT_4 yCalCorrection = 0;
	INT_4 bCalCorrection = 0;
	INT_4 nOxMetCount = 0;
	INT_4 cOxMetCount = 0;
	INT_4 bSingleAACount = 0, ySingleAACount = 0, bSingleAASeries = 0, ySingleAASeries = 0;
	char test, twoAAExtension, maxCharge;
	BOOLEAN monoToAvYSwitch = TRUE;	/*Used to recalculate for average masses.*/
	BOOLEAN avToMonoBSwitch = FALSE;	/*Used to recalculate for average masses.*/
	BOOLEAN twoAANTerm;
	char testForPro;
	REAL_4 currentIonFound;
			
/*	Initialize a few variables.*/
	bCount        = 0;	/*Current number of b ions in a row.*/
	bCountStringent = 0;	/*A more stringent count used for quality assessment of the spectrum*/
	yCount        = 0;	/*Ditto.*/
	yCountStringent = 0;
	bSeries       = 0;	/*The greatest number of b ions in a row for a given sequence.*/
	bSeriesStringent = 0;
	ySeries       = 0;	/*The greatest number of y ions in a row for a given sequence.*/
	ySeriesStringent = 0;
	skipOneY	  = 1;	/*Allows for one missed y ion in a series w/o resetting to zero.*/
	skipOneB	  = 1;	/*Allows for one missed b ion in a series w/o resetting to zero.*/
	cleavageSites = 0;	/*Used to count the number of times a y or b ion of any charge state 
						delineates a sequence. */
	gCleavageSiteStringent = 0;
	precursor 	  = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) 
				    / gParam.chargeState;
				    
	lossOf64 = gElementMass_x100[CARBON] + gElementMass_x100[SULFUR] + 
				gElementMass_x100[OXYGEN] + 4 * gElementMass_x100[HYDROGEN];
				    
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}

/*	
	Initialize the starting b ion mass (acetylated, etc).
*/

	bCalStart = gParam.modifiedNTerm + 0.5;
	/*Determine the correction factor for high accuracy*/
	i = gParam.modifiedNTerm * 10 + 0.5;
	bCalCorrection = i - bCalStart * 10;
	
	yCalStart = gParam.modifiedCTerm + + (2 * gElementMass_x100[HYDROGEN]) + 0.5;
	i = (gParam.modifiedCTerm + (gMultiplier * gElementMass[HYDROGEN] * 2)) * 10 + 0.5;
	yCalCorrection = i - yCalStart * 10;
	
/*	
	Count the number of charged residues.  nChargeCount is one more than the number of 
	charged residues found in an N-terminal fragment.  cChargeCount is one more than the 
	number of charged residues found in a C-terminal fragment.
*/
	nChargeCount = FindNCharge(sequence, seqLength);
	cChargeCount = 1;
	
/*
	Determine if oxidized methionine is present.  nOxMetCount and cOxMetCount are integer values,
	depending on if the b and y ions contain an oxidized methionine.
*/

	nOxMetCount = FindNOxMet(sequence, seqLength);
	cOxMetCount = 0;
	
/*	
	Figure out if the N-terminus is a two amino acid extension (TRUE or FALSE).
*/
	twoAANTerm = TwoAAExtFinder(sequence, 0);

	for(i = (seqLength - 1); i > 0; i--)	/*Don't do this loop for i = 0 (doesnt make sense).*/
	{
		
/*	
	Initialize some variables for this 'for' loop.
*/
		bOrYIon = 0;	/*If any number of b or y ions of any charge are found, 
						then this equals one.  Otherwise, it stays at zero.*/
		bIon = 0;	/*If a b ion is found this equals one.*/
		yIon = 0;	/*If a y ion is found this equals one.*/
		
/*	
	Figure out if this is a two amino acid extension (TRUE or FALSE).
*/
		twoAAExtension = TwoAAExtFinder(sequence, i);

/*
	Calculate the singly charged y ion mass.
*/
		yCal = YCalculator(i, sequence, seqLength, yCalStart, yCalCorrection);	
			
/*
	Calculate the singly charged b ion mass.
*/
		bCal = BCalculator(i, sequence, bCalStart, bCalCorrection);
				
/*
	If the N-terminus is not a two amino acid extension, its unlikely that one can find
	a b1 ion, so give bCal an irrelevant value, so that an accidental match to b1 is not
	made.
*/
		if(i == 1)
		{
			/*find out if first amino acid is a twoAA extension*/
			twoAAExtension = TwoAAExtFinder(sequence, i - 1);
			if(!twoAAExtension)
			{
				bCal = 10;
			}
			/*get the twoAAextension for the second position again*/
			twoAAExtension = TwoAAExtFinder(sequence, i);
		}
		
/*
	Readjust the number of oxidized methionines.
*/

		if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0)	/*mass accuracy sufficient
																	to determine oxMet*/
		{
			if(sequence[i] >= gMonoMass_x100[9] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[9] + gToleranceNarrow)
			{
				nOxMetCount--;	/*b ions lost a oxMet*/
				cOxMetCount++;	/*y ions gained a oxMet*/
			}
			if(nOxMetCount < 0 || cOxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}
		else	/*mass accuracy not sufficient to differentiate oxMet from Phe*/
		{
			if(sequence[i] >= gMonoMass_x100[F] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[F] + gToleranceNarrow)
			{
				nOxMetCount--;	/*b ions lost a oxMet*/
				cOxMetCount++;	/*y ions gained a oxMet*/
			}
			if(nOxMetCount < 0 || cOxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}

/*	
	Readjust the number of charges in the C- and N-terminii.
*/
		if((sequence[i] >= gMonoMass_x100[R] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[R] + gToleranceWide) || 
			(sequence[i] >= gMonoMass_x100[H] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[H] + gToleranceWide) ||
			(sequence[i] >= gMonoMass_x100[K] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[K] + gToleranceWide))
		{
			cChargeCount += 1;
			nChargeCount -= 1;
		}
		else	/*Check to see if its a two amino acid combo that could contain Arg, His, or Lys.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] >= gArgPlus[j] - gToleranceWide 
					&& sequence[i] <= gArgPlus[j] + gToleranceWide) || 
					(sequence[i] >= gHisPlus[j] - gToleranceWide 
					&& sequence[i] <= gHisPlus[j] + gToleranceWide) ||
					(sequence[i] >= gLysPlus[j] - gToleranceWide 
					&& sequence[i] <= gLysPlus[j] + gToleranceWide))
				{
					cChargeCount += 1;
					nChargeCount -= 1;
					break;
				}
			}
		}
		
/*
	Check each charge state up to the parent ion charge.
*/
		for(j = 1; j <= maxCharge; j++) 
		{

/*	
	Initialize variables within this j loop.
*/
			test = FALSE;	/*Used to test if b, a, or y ions are found before looking 
							for the corresponding losses of ammonia or water.*/
			aIonMass = (bCal - gCO + (j * gElementMass_x100[HYDROGEN]) 
						- gElementMass_x100[HYDROGEN]) / j;
				aIonMassMinErr = aIonMass - gToleranceWide;
				aIonMassPlusErr = aIonMass + gToleranceWide;
			bIonMass = (bCal + (j * gElementMass_x100[HYDROGEN]) - 
						gElementMass_x100[HYDROGEN]) / j;
				bIonMassMinErr = bIonMass - gToleranceWide;
				bIonMassPlusErr = bIonMass + gToleranceWide;
			yIonMass = (yCal + (j * gElementMass_x100[HYDROGEN]) 
						- gElementMass_x100[HYDROGEN]) / j;
				yIonMassMinErr = yIonMass - gToleranceWide;
				yIonMassPlusErr = yIonMass + gToleranceWide;
			
/*	Search for b ions.*/
			if((bIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				k = fragNum - 1;
				while(fragMOverZ[k] >= bIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= bIonMassPlusErr)
					{
						if(nChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							bOrYIon = 1;	/*A b or y ion has been identified.*/
							test = TRUE;	/*A b ion of charge j has been identified.*/
							bIon = 1;		/*A b ion is present.*/
							massDiff = abs(bIonMass - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							
	
/*	
	If the charge under consideration "j" is equal to the precursor charge, and if that precursor
	charge is greater than one, then the bion intensity is attenuated.  Alternatively, if
	the b ion mass is greater than the precusor, while at the same time the number of basic
	residues in the b ion is less than the number of charges on the precursor, then that is also
	grounds for reducing the influence of the ion.
*/
							if((j == gParam.chargeState && gParam.chargeState > 1) || 
								(bIonMass > precursor && nChargeCount < gParam.chargeState))
							{
								if(gParam.fragmentPattern != 'L')
								{
									ionFound[k] = ionFound[k] * HIGH_MASS_B_ION_MULTIPLIER;
								}
							}
							if(twoAAExtension)
							{
								ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
							}
							if(currentIonFound > ionFound[k])	/*if old number is bigger than new*/
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								ionType[k] = 2;	/*b ions are type 2*/
							}
							if(j < gParam.chargeState && (bIonMass < precursor 
								|| gParam.fragmentPattern == 'L'))
							{
								bFound[k] = 1;
								byError[k] = AssignError(byError[k], bIonMass, fragMOverZ[k]);
							}
						}
					}
					k--;
				}

/*	
*	Search for b minus ammonia or water.  The index value of k is carried over from the 
*	b ion search, since
*	the while loop breaks out once the mass value is less than the calculated b value 
*	(ie, k is close to 
*	b minus ammonia or water.
*/
				if(!test)
				{
					if(sequence[0] <= gMonoMass_x100[K] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[K] - gParam.fragmentErr
						&& gParam.fragmentErr > 0.04 * gMultiplier)
					{
						test = TRUE;
					}
					if(sequence[0] <= gMonoMass_x100[Q] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[Q] - gParam.fragmentErr)
					{
						test = TRUE;
					}
					if(!test)
					{
						for(m = 0; m < gAminoAcidNumber; m++)
						{
							if(sequence[0] <= gGlnPlus[m] + gParam.fragmentErr &&
								sequence[0] >= gGlnPlus[m] - gParam.fragmentErr)
							{
								test = TRUE;
								break;
							}
							if(sequence[0] + sequence[1] <= gGlnPlus[m] + gParam.fragmentErr &&
								sequence[0] + sequence[1] >= gGlnPlus[m] - gParam.fragmentErr)
							{
								test = TRUE;
								break;
							}
						}
					}
				}
				
				if(!test)
				{
					if(sequence[0] <= gMonoMass_x100[E] + gParam.fragmentErr && 
						sequence[0] >= gMonoMass_x100[E] - gParam.fragmentErr)
					{
						test = TRUE;
					}
					
					if(!test)
					{
						for(m = 0; m < gAminoAcidNumber; m++)
						{
							if(sequence[0] <= gGluPlus[m] + gParam.fragmentErr &&
								sequence[0] >= gGluPlus[m] - gParam.fragmentErr)
							{
								test = TRUE;
								break;
							}
							if(sequence[0] + sequence[1] <= gGluPlus[m] + gParam.fragmentErr &&
								sequence[0] + sequence[1] >= gGluPlus[m] - gParam.fragmentErr)
							{
								test = TRUE;
								break;
							}
						}
					}
				}
					
				if(test || argPresent)
				{
				
/*	Calculate the b minus water and b minus ammonia values.*/
					bMinW = (bCal - gWater + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						bMinWMinErr = bMinW - gToleranceWide;
					bMinA = (bCal - gAmmonia + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						bMinAPlusErr = bMinA + gToleranceWide;
					
					if(nChargeCount >= j)
					{
						while(fragMOverZ[k] >= bMinWMinErr && k >= 0)
						{
							if(fragMOverZ[k] <= bMinAPlusErr)
							{
								massDiffW = abs(bMinW - fragMOverZ[k]);
								massDiffA = abs(bMinA - fragMOverZ[k]);
								if(massDiffW < massDiffA)
								{
									massDiff = massDiffW;
								}
								else
								{
									massDiff = massDiffA;
								}
								currentIonFound = ionFound[k];
								ionFound[k] = CalcIonFound(ionFound[k], massDiff);
								if(gParam.chargeState > 1)
								{
									ionFound[k] = ionFound[k] * NEUTRAL_LOSS_MULTIPLIER;
									if((j == gParam.chargeState && gParam.chargeState > 1) || 
										(bMinWMinErr > precursor && nChargeCount < gParam.chargeState))
									{
										if(gParam.fragmentPattern != 'L')
										{
											ionFound[k] = ionFound[k] * HIGH_MASS_B_ION_MULTIPLIER;
										}
									}
								}
								if(twoAAExtension)
								{
									ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
								}
								if(currentIonFound > ionFound[k])
								{
									ionFound[k] = currentIonFound;
								}
								else
								{
									ionType[k] = 3;	/*b-17/18 ions are type 2*/
								}
	
							}
							k--;
						}
					}
				}
			}	/*if((bIonMass * j) > ((j-1) * 50000))*/
				
/*	
*	Search for the a ions.  
*	The k index value is retained from the b ion search, and the search for a continues 
*	downward.
*/
			/*Make sure there is enough mass to hold the charge, and that there is a b ion.*/
			if((aIonMass * j) > ((j-1) * 400 * gMultiplier) && test) 
			{
				test = FALSE;
				while(fragMOverZ[k] >= aIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= aIonMassPlusErr)
					{
						if(nChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							test = TRUE;	/*An a ion of charge j has been identified.*/
							massDiff = abs(aIonMass - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							if((i > 2 && !twoAANTerm) || (i > 1 && twoAANTerm))
							{
								if((j == gParam.chargeState && gParam.chargeState > 1) || 
									(aIonMass > precursor && nChargeCount < gParam.chargeState))
								{
									ionFound[k] = ionFound[k] * HIGH_MASS_A_ION_MULTIPLIER
																* HIGH_MASS_A_ION_MULTIPLIER;
								}
								else
								{
									ionFound[k] = ionFound[k] * HIGH_MASS_A_ION_MULTIPLIER;
								}
							}
							if(twoAAExtension)
							{
								ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
							}
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								ionType[k] = 4;	/*a ions are type 2*/
							}

						}
					}
					k--;
				}

/*	
*	Search for a minus ammonia or water.  The index value of k is carried over from the 
*	a ion search, since
*	the while loop breaks out once the mass value is less than the calculated a value 
*	(ie, k is now close to 
*	a minus ammonia or water.
*/
				if(test || argPresent)
				{
				
/*	Calculate the a minus water and a minus ammonia values.*/
					aMinW = (bCal - gCO - gWater + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						aMinWMinErr = aMinW - gToleranceWide;
					aMinA = (bCal - gCO - gAmmonia + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						aMinAPlusErr = aMinA + gToleranceWide;
					if(nChargeCount >= j)
					{
						while(fragMOverZ[k] >= aMinWMinErr && k >= 0)
						{
							if(fragMOverZ[k] <= aMinAPlusErr)
							{
								massDiffW = abs(aMinW - fragMOverZ[k]);
								massDiffA = abs(aMinA - fragMOverZ[k]);
								if(massDiffW < massDiffA)
								{
									massDiff = massDiffW;
								}
								else
								{
									massDiff = massDiffA;
								}
								currentIonFound = ionFound[k];
								ionFound[k] = CalcIonFound(ionFound[k], massDiff);
								if((i > 2 && !twoAANTerm) || (i > 1 && twoAANTerm) || j != 1)
								{
									ionFound[k] = ionFound[k] * HIGH_MASS_A_ION_MULTIPLIER;
								}
								if(gParam.chargeState > 1)
								{
									ionFound[k] = ionFound[k] * NEUTRAL_LOSS_MULTIPLIER;
								}
								if(twoAAExtension)
								{
									ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
								}
								if(currentIonFound > ionFound[k])
								{
									ionFound[k] = currentIonFound;
								}
								else
								{
									ionType[k] = 5;	/*a-17/18 ions are type 2*/
								}
							}
							k--;
						}
					}
				}
			}	/*if((aIonMass * j) > ((j-1) * 50000))*/
			
/*	
*	Search for the b minus 64 ions if oxMet is present.  
*	The k index value is retained from the b and a ion search, and continues 
*	downward.
*/
			if(test && nOxMetCount) /*there's a b ion and oxMet is in the b ion*/
			{
				bMin64 = (bCal - lossOf64 + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
				bMin64MinErr = bMin64 - gToleranceWide;
				bMin64PlusErr = bMin64 + gToleranceWide;
					
				while(fragMOverZ[k] >= bMin64MinErr && k >= 0)
				{
					if(fragMOverZ[k] <= bMin64PlusErr)
					{
						massDiff = abs(bMin64 - fragMOverZ[k]);
						currentIonFound = ionFound[k];
						ionFound[k] = CalcIonFound(ionFound[k], massDiff);
						if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0)
						{
							ionFound[k] = ionFound[k] * OXMET_MULTIPLIER;
						}
						else
						{
							ionFound[k] = ionFound[k] * PHE_MULTIPLIER;
						}
						if(currentIonFound > ionFound[k])
						{
							ionFound[k] = currentIonFound;
						}
						else
						{
							ionType[k] = 6;	/*losses from oxidized Met ions are type 2*/
						}
					}
					k--;
				}
			}
/*	Search for the y ion values.*/
			if((yIonMass * j) > ((j-1) * 400 * gMultiplier)) /*Make sure there is enough mass to hold 
													the charge.*/
			{
				test = FALSE;
				k = fragNum - 1;
				while(fragMOverZ[k] >= yIonMassMinErr && k >= 0)
				{
					if(fragMOverZ[k] <= yIonMassPlusErr)
					{
						if(cChargeCount >= j)	/*Make sure enough charges can be attached.*/
						{
							bOrYIon = 1;	/*A b or y ion has been identified.*/
							test = TRUE;	/*A y ion of charge j has been identified.*/
							yIon = 1;		/*A y ion is present.*/
							massDiff = abs(yIonMass - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							if(j == gParam.chargeState && gParam.chargeState > 1)
							{
								if((((i > 2 && !twoAANTerm) || (i > 1 && twoAANTerm)) 
										&& gParam.fragmentPattern != 'L') || 
										(i > (INT_4)seqLength / 3 && gParam.fragmentPattern == 'L'))
								
								{
									ionFound[k] = ionFound[k] * HIGH_CHARGE_Y_ION_MULTIPLIER;
								}
							}
							if(twoAAExtension)
							{	/*see if the twoaa ext could contain pro*/
								testForPro = FALSE;
								for(m = 0; m < gAminoAcidNumber; m++)
								{
									if(sequence[i] <= gProPlus[m] + gToleranceWide &&
										sequence[i] >= gProPlus[m] - gToleranceWide)
									{
										testForPro = TRUE;
									}
								}
								if(!testForPro)	/*if could contain pro then don't attenuate*/
								{
									ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
								}
							}
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								ionType[k] = 7;	/*y ions are type 7*/
							}
							if(j < gParam.chargeState)
							{
								yFound[k] = 1;
								byError[k] = AssignError(byError[k], yIonMass, fragMOverZ[k]);
							}
						}
					}
					k--;
				}

/*	
*	Search for y minus ammonia or water.  The index value of k is carried over from the 
*	y ion search, since
*	the while loop breaks out once the mass value is less than the calculated y value 
*	(ie, k is close to 
*	y minus ammonia or water.
*/
				if(test || argPresent)
				{
				
/*	Calculate the y minus water and y minus ammonia values.*/
					yMinW = (yCal - gWater + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						yMinWMinErr = yMinW - gToleranceWide;
					yMinA = (yCal - gAmmonia + (j * gElementMass_x100[HYDROGEN]) 
							- gElementMass_x100[HYDROGEN]) / j;
						yMinAPlusErr = yMinA + gToleranceWide;
					if(cChargeCount >= j)
					{
						while(fragMOverZ[k] >= yMinWMinErr && k >= 0)
						{
							if(fragMOverZ[k] <= yMinAPlusErr)
							{
								massDiffW = abs(yMinW - fragMOverZ[k]);
								massDiffA = abs(yMinA - fragMOverZ[k]);
								if(massDiffW < massDiffA)
								{
									massDiff = massDiffW;
								}
								else
								{
									massDiff = massDiffA;
								}
								currentIonFound = ionFound[k];
								ionFound[k] = CalcIonFound(ionFound[k], massDiff);
								if(gParam.chargeState > 1)
								{
									ionFound[k] = ionFound[k] * NEUTRAL_LOSS_MULTIPLIER;
									if(j == gParam.chargeState)
									{
										ionFound[k] = ionFound[k] * HIGH_CHARGE_Y_ION_MULTIPLIER;
									}
								}
								if(twoAAExtension)
								{
									ionFound[k] = ionFound[k] * TWO_AA_EXTENSION_MULTIPLIER;
								}
								if(currentIonFound > ionFound[k])
								{
									ionFound[k] = currentIonFound;
								}
								else
								{
									ionType[k] = 8;	/*y-17/18 ions are type 8*/
								}
							}
							k--;
						}
					}
					
/*	
*	Search for the b minus 64 ions if oxMet is present.  
*	The k index value is retained from the b and a ion search, and continues 
*	downward.
*/
					if(test && cOxMetCount) /*there's a b ion and oxMet is in the b ion*/
					{
						yMin64 = (yCal - lossOf64 + (j * gElementMass_x100[HYDROGEN]) 
									- gElementMass_x100[HYDROGEN]) / j;
						yMin64MinErr = yMin64 - gToleranceWide;
						yMin64PlusErr = yMin64 + gToleranceWide;
					
						while(fragMOverZ[k] >= yMin64MinErr && k >= 0)
						{
							if(fragMOverZ[k] <= yMin64PlusErr)
							{
								massDiff = abs(yMin64 - fragMOverZ[k]);
								currentIonFound = ionFound[k];
								ionFound[k] = CalcIonFound(ionFound[k], massDiff);
								if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0)
								{
									ionFound[k] = ionFound[k] * OXMET_MULTIPLIER;
								}
								else
								{
									ionFound[k] = ionFound[k] * PHE_MULTIPLIER;
								}
								if(currentIonFound > ionFound[k])
								{
									ionFound[k] = currentIonFound;
								}
								else
								{
									ionType[k] = 9;	/*losses from oxMet ions are type 2*/
								}
							}
							k--;
						}
					}
					
				}
			}	/*if((yIonMass * j) > ((j-1) * 50000))*/
		}	/*for j*/
		if(bIon)
		{
			bCount++;	/*If there was a b ion, then increment by one.*/
			skipOneB = 1;
			bCountStringent++;	/*this count doesn't allow for gaps in the series*/
		}
		else
		{
			if(skipOneB == 0)
			{
				bCount = 0;	/*Otherwise, reset the counting of b ions to zero.*/
			}
			skipOneB = 0;	/*Allow one missing b ion before resetting bCount*/
			bCountStringent = 0;
		}
		/*Sometimes a K or R is forced into the C-terminus even if no fragments are found.
		If the C-terminus is not a gap, then increment bCount, but don't reset the skipOneB,
		since if the penultimate amino acid lacks cleavage then it should be trashed.*/
		if(twoAAExtension == FALSE && i == seqLength - 1 && !bIon)
		{
			bCount++;
		}
		/*bSingleAACount counts the contiguous b ions differing by single amino acids.*/
		if(bIon && !twoAAExtension)
		{
			bSingleAACount++;
		}
		else if(bIon && twoAAExtension)
		{
			bSingleAACount = 1;	/*reset to one so that C-terminal cleavage next to dipeptide is counted*/
		}
		else
		{
			bSingleAACount = 0;
		}
		
		if(yIon)
		{
			yCount++;
			skipOneY = 1;
			yCountStringent++;	/*this count doesn't allow for gaps in the series*/
		}
		else
		{
			if(skipOneY == 0)
			{
				yCount = 0;
			}
			skipOneY = 0;
			yCountStringent = 0;
		}
		if(twoAAExtension == FALSE && i == seqLength - 1 && !yIon)
		{
			yCount++;
		}
		/*ySingleAACount counts the contiguous y ions differing by single amino acids.*/
		if(yIon && !twoAAExtension)	/*if there's a y ion and its a single amino acid*/
		{
			ySingleAACount++;
		}
		else if(yIon && twoAAExtension)	/*if there's a y ion and its a two aa extension*/
		{
			ySingleAACount = 1;	/*reset to one so that C-terminal cleavage next to dipeptide is counted*/
		}
		else
		{
			ySingleAACount = 0;
		}
		
		if(bCount > bSeries)
		{
			bSeries = bCount;	/*Don't forget what the longest continuous b series was.*/
		}
		if(yCount > ySeries)
		{
			ySeries = yCount;	/*Don't forget what the longest continuous y series was.*/
		}
		if(bCountStringent > bSeriesStringent)
		{
			bSeriesStringent = bCountStringent;
		}
		if(yCountStringent > ySeriesStringent)
		{
			ySeriesStringent = yCountStringent;
		}
		/*bSingleAASeries and ySingleAASeries keep track of the longest stretch of single aa seq*/
		if(bSingleAACount > bSingleAASeries)
		{
			bSingleAASeries = bSingleAACount;
		}
		if(ySingleAACount > ySingleAASeries)
		{
			ySingleAASeries = ySingleAACount;
		}
		/*if(gParam.fragmentPattern == 'L' && (i < 3 || i > (seqLength - 3)))
		{
			cleavageSites++;*/	/*For LCQ data, the ends are not well established and should
								not be counted against sequences*/
		/*}
		else
		{
			cleavageSites += bOrYIon;*/	/*Count the number of times a b or y ion define a 
									cleavage site.*/
		/*}*/

		
	}	/*for i*/
	
/*
*	For TSQ data I found that I get better scoring if I look for contiguous series of either
*	b or y ions.  For LCQ data, because of the missing low mass end, I cannot expect a 
*	contiguous series, especially since the high mass b ions will compensate for the loss of 
*	low mass y ions.
*/

	if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q' 
		|| gParam.fragmentPattern == 'L')
	{
		if(ySeries > bSeries)
		{
			cleavageSites = ySeries;
			gCleavageSiteStringent = ySeriesStringent;
		}
		else
		{
			cleavageSites = bSeries;
			gCleavageSiteStringent = bSeriesStringent;
		}
		/*This is the number of cleavage sites defining a contiguous series of single aa's*/
		if(ySingleAASeries > bSingleAASeries)
		{
			gSingleAACleavageSites = ySingleAASeries;
		}
		else
		{
			gSingleAACleavageSites = bSingleAASeries;
		}
	}
		
	return(cleavageSites);
}

/******************************InternalFrag************************************************
*
*	InternalFrag identifies internal fragment ions (where neither the C- nor N-terminii 
*	are present).  Included are the losses of CO, water, and ammonia from the usual b-type 
*	internal fragment.  The input is as described in the documentation for the function 
*	PEFragments, and it returns nothing.  Only the array ionFound is modified.  Only 
*	singly-charged ions are considered here.
*/
void InternalFrag(REAL_4 *ionFound, INT_4 *fragMOverZ, 
				INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *ionType)
{
	INT_4 i, j, k, intFrag, intFragMinErr, intFragPlusErr, saveIndex;
	INT_4 intFragMinW, intFragMinWMinErr;
	INT_4 intFragMinA, intFragMinAPlusErr;
	INT_4 intFragMinCO, intFragMinCOMinErr, intFragMinCOPlusErr;
	INT_4 massDiff, massDiffW, massDiffA, precursor;
	char intFragFound;
	REAL_4 currentIonFound;
	
	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) 
				/ gParam.chargeState;

	if(seqLength >= 4)	/*Sequences less than 4 amino acids cannot have internal fragment ions.*/
	{
		for(i = 1; i < (seqLength - 3); i++)	/*This is the N-terminus of the fragment.*/
		{
			for(j = (i + 1); j < (seqLength - 2); j++)	/*C-terminus of the fragment.*/
			{
				intFragFound = 0;	/*This tests to see if the standard int frag is present; 
										if it is present, then the losses of CO, water, 
										and ammonia are also searched for.*/
										
				intFrag = gElementMass_x100[HYDROGEN];	/*Calc the mass of the int frag.*/
				for(k = i; k <= j; k++)
				{
					intFrag += sequence[k];
				}
				intFragMinErr = intFrag - gToleranceWide;
				intFragPlusErr = intFrag + gToleranceWide;
				
				if(intFragMinErr < precursor)	/*Only count those matches where the internal frag
												mass is less than the precursor m/z value.*/
				{
					k = 0;
					while(fragMOverZ[k] <= intFragPlusErr && k < fragNum)	/*Look for this ion.*/
					{
						if(fragMOverZ[k] >= intFragMinErr)
						{
							currentIonFound = ionFound[k];
							massDiff = abs(intFrag - fragMOverZ[k]);
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							/*Attenuate the internal ion intensity, unless its a short one w/ P at N-term*/
							if(sequence[i] != gGapList[P] || j > i + 2)
							{
								ionFound[k] = ionFound[k] * INTERNAL_FRAG_MULTIPLIER ;
							}
							intFragFound = 1;	/*Its been found.*/
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								if(sequence[i] != gGapList[P] || j > i + 2)
								{
									ionType[k] = 10;	/*internal frag w/o Pro*/
								}
								else
								{
									ionType[k] = 11;	/*internal frag w/ Pro*/
								}
							}
						}
						k++;
					}
					saveIndex = k;	/*Start looking for the next fragment types using this index.*/
				}
				if(intFragFound)
				{
/*	Calculate the various losses and plus/minus tolerances.	*/
					intFragMinW = intFrag - gWater;
						intFragMinWMinErr = intFragMinW - gToleranceWide;
					intFragMinA = intFrag - gAmmonia;
						intFragMinAPlusErr = intFragMinA + gToleranceWide;
					intFragMinCO = intFrag - gCO;
						intFragMinCOMinErr = intFragMinCO - gToleranceWide;
						intFragMinCOPlusErr = intFragMinCO + gToleranceWide;
					
/*	Calculate internal fragments minus water or ammonia.	*/
					k = saveIndex;
					while(fragMOverZ[k] >= intFragMinWMinErr && k >= 0)
					{
						if(fragMOverZ[k] <= intFragMinAPlusErr)
						{
							massDiffW = abs(intFragMinW - fragMOverZ[k]);
							massDiffA = abs(intFragMinA - fragMOverZ[k]);
							if(massDiffW < massDiffA)
							{
								massDiff = massDiffW;
							}
							else
							{
								massDiff = massDiffA;
							}
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							ionFound[k] = ionFound[k] * INTERNAL_FRAG_MULTIPLIER 
														* NEUTRAL_LOSS_MULTIPLIER;
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								ionType[k] = 12;	/*internal frag minus 17 or 18*/
							}

						}
						k--;
					}
					
/*	Calculate internal fragments minus carbon monoxide.*/
					k = saveIndex;
					while(fragMOverZ[k] >= intFragMinCOMinErr && k >= 0)
					{
						if(fragMOverZ[k] <= intFragMinCOPlusErr)
						{
							massDiff = abs(intFragMinCO - fragMOverZ[k]);
							currentIonFound = ionFound[k];
							ionFound[k] = CalcIonFound(ionFound[k], massDiff);
							ionFound[k] = ionFound[k] * INTERNAL_FRAG_MULTIPLIER 
														* NEUTRAL_LOSS_MULTIPLIER;
							if(currentIonFound > ionFound[k])
							{
								ionFound[k] = currentIonFound;
							}
							else
							{
								ionType[k] = 13;	/*internal frag minus CO*/
							}
						}
						k--;
					}
				}	/*if intFragFound*/
			}	/*for j*/
		}	/*for i*/
	}	/*if sequence is INT_4 enough*/

	return;	
}

/*******************************ArgIons************************************
*
*	ArgIons inputs the usual stuff for identifying ions (see PEFragments documentation).  
*	It returns a char value of 0 or 1, depending whether Arg is present in the current 
*	sequence.  This function looks for an ion I call b - OH, which is unique to singly-charged 
*	precursor ions containing arginine.  The ion b-OH-NH3 is also identified.  Only 
*	singly-charged fragments are considered here.
*/
char ArgIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
			 INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, j, bOH, bOHMinErr, bOHPlusErr, bOHMinAmm, bOHMinAmmMinErr, 
				bOHMinAmmPlusErr, massDiff, argCount;
	char argPresent = FALSE;
	char bOHPresent = 0;
	
	argCount = 0;	/*counts the number of arginines in the sequence*/
	
/*	Determine if Arg is in the sequence as a single amino acid. */
	for(i = 0; i < seqLength; i++)
	{
		if((sequence[i] <= gMonoMass_x100[R] + gToleranceWide) &&
			(sequence[i] >= gMonoMass_x100[R] -gToleranceWide))
		{
	/*		argPresent = TRUE;*/
			argCount++;	
		}
	}
	
	if(argPresent == FALSE)	/*Check to see if Arg is present as a two amino acid gap.*/
	{
		for(i = 0; i < seqLength; i++)	/*For each "amino acid" in the array sequence.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] <= gArgPlus[j] + gToleranceWide) &&
					(sequence[i] >= gArgPlus[j] - gToleranceWide))	/*Is the amino acid one of the two 
															amino acid combinations that contain 
															arginine?*/
				{
		/*			argPresent = TRUE;*/
					argCount++;
				/*	break;*/
				}
			}
		/*	if(argPresent)
		//	{
		//		break;*/	/*If I found arginine, why look for more?*/
		/*	}*/
		}
	}

	/*If the number of arginines exceeds the number of protons, then its a non-mobile proton situation, which
	is designated by argPresent becoming TRUE.  Also, the bOH ions are counted here*/
	if(argCount >= gParam.chargeState)
	{
		argPresent = TRUE;		
		bOH = gParam.peptideMW + gElementMass_x100[HYDROGEN] - sequence[seqLength - 1];
		bOHMinErr = bOH - gToleranceWide;
		bOHPlusErr = bOH + gToleranceWide;
		
		i = fragNum - 1;
		while(fragMOverZ[i] >= bOHMinErr && i >= 0)
		{
			if(fragMOverZ[i] <= bOHPlusErr)
			{
				massDiff = abs(bOH - fragMOverZ[i]);
				ionFound[i] = CalcIonFound(ionFound[i], massDiff);
				bOHPresent = 1;
			}
			i--;
		}
		if(bOHPresent)
		{
			bOHMinAmm = bOH - gAmmonia;
			bOHMinAmmMinErr = bOHMinAmm - gToleranceWide;
			bOHMinAmmPlusErr = bOHMinAmm + gToleranceWide;
			
			i = fragNum - 1;
			while(fragMOverZ[i] >= bOHMinAmmMinErr && i >= 0)
			{
				if(fragMOverZ[i] <= bOHMinAmmPlusErr)
				{
					massDiff = abs(bOHMinAmm - fragMOverZ[i]);
					ionFound[i] = CalcIonFound(ionFound[i], massDiff);
				}
				i--;
			}
		}
	}

	return(argPresent);
}							
									
/****************************LoadSequence***********************************************
*
*	LoadSequence inputs a INT_4 array called "sequence", a pointer to a INT_4 
*	"seqLength", and a pointer to a struct of type SequenceData.  It takes the 
*	INT_4 sequence found in the "peptide" field of the struct, and puts a sequence 
*	of integers in the "sequence" array.  This sequence of integers corresponds to the monoisotopic
*	residue mass of single amino acids or pairs of amino acids.
*	number for the amino acids found in gSingAA.
*/

void LoadSequence(INT_4 *sequence, INT_4 *seqLength, struct Sequence *currSeqPtr)
{
	INT_4 i, j, testMass;
	char test;
	
	*seqLength = currSeqPtr->peptideLength;
	
	for(i = 0; i < *seqLength; i++)
	{
		test = TRUE;
		/*make sure the correct mass is used*/
		testMass = currSeqPtr->peptide[i];
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			/*if((testMass <= gGapList[j] + gToleranceNarrow)
				&& (testMass >= gGapList[j] - gToleranceNarrow))*/
			if(testMass == gGapList[j])
			{
				sequence[i] = gMonoMass_x100[j];
				test = FALSE;
				break;
			}
		}
		if(test)	/*if two aa extension*/
		{
			sequence[i] = testMass;
		}
	}
	
		
	return;
}

/********************************WaterLoss*************************************************
*
*	WaterLoss identifies fragment ions that are due to losses of two waters, two ammonias 
*	or one of each.  The input is the ionFound array (0 = not identified, 1 = identified), 
*	fragNum (fragment ion number), fragmentErr (fragment ion tolerance), fragMOverZ array 
*	(of fragment ion masses), and chargeState (the precursor ion charge state).  It 
*	returns nothing, but it does modify the ionFound array if ions match w/ calculated values.
*/

void WaterLoss(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, INT_4 *fragIntensity)
{
	INT_4 precurMin2W, precurMinWA, precurMin2A, precurMinW, precurMinA;
	INT_4 precurMin2WMinErr, precurMin2WPlusErr, precurMinWAMinErr, 
				precurMinWAPlusErr, precurMin2AMinErr, precurMin2APlusErr,
				precurMinWPlusErr, precurMinWMinErr,
				precurMinAPlusErr, precurMinAMinErr;
	INT_4 i;

	precurMin2W = (gParam.peptideMW - gWater - gWater + 
	              (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	precurMinWA = (gParam.peptideMW - gWater - gAmmonia +
	              (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	precurMin2A = (gParam.peptideMW - gAmmonia - gAmmonia + 
	              (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	precurMinW = (gParam.peptideMW - gWater + 
	             (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	precurMinA = (gParam.peptideMW - gAmmonia + 
	             (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;

	precurMin2WMinErr = precurMin2W - gToleranceWide;
	precurMin2WPlusErr = precurMin2W + gToleranceWide;
	precurMinWAMinErr = precurMinWA - gToleranceWide;
	precurMinWAPlusErr = precurMinWA + gToleranceWide;
	precurMin2AMinErr = precurMin2A - gToleranceWide;
	precurMin2APlusErr = precurMin2A + gToleranceWide;
	precurMinWPlusErr = precurMinW + gToleranceWide;
	precurMinWMinErr = precurMinW - gToleranceWide;
	precurMinAPlusErr = precurMinA + gToleranceWide;
	precurMinAMinErr = precurMinA - gToleranceWide;
	
	
	for(i = 0; i < fragNum; i++)
	{
		if(fragMOverZ[i] <= precurMin2APlusErr)
		{
			if(fragMOverZ[i] >= precurMin2WMinErr)
			{
				if(((fragMOverZ[i] <= precurMin2WPlusErr) && (fragMOverZ[i] >= precurMin2WMinErr))
				|| ((fragMOverZ[i] <= precurMinWAPlusErr) && (fragMOverZ[i] >= precurMinWAMinErr))
				|| ((fragMOverZ[i] <= precurMin2APlusErr) && (fragMOverZ[i] >= precurMin2AMinErr))
				|| ((fragMOverZ[i] <= precurMinWPlusErr) && (fragMOverZ[i] >= precurMinWMinErr))
				|| ((fragMOverZ[i] <= precurMinAPlusErr) && (fragMOverZ[i] >= precurMinAMinErr)))
				{
					ionFound[i] = 1;
				}
			}
		}
	}
/*	If an ion was converted to intensity of zero in the function TotalIntensity, it is 
	identified as found here.*/
	for(i = 0; i < fragNum; i++)
	{
		if(fragIntensity[i] == 0)
		{
			ionFound[i] = 1;
		}
	}
	return;
}

/****************************ScoreC1*************************************************
*
*	If the first position is a two amino acid residue, then check to see if it could
*	be Gln and if the remaining mass corresponds to an amino acid.  If so, define the
*	residue as two separate ones X and Q.  If the N-terminus is a single amino acid,then
*	check if the next one is Q.  Score the c1 ion.
*/
INT_4 ScoreC1(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, 
				 INT_4 *sequence, INT_4 seqLength)
{
	INT_4	i, c1Ion, testMass, massDiff;
	REAL_4	currentIonFound;
	BOOLEAN	oneAA;
	
/*	Initalize*/
	oneAA				= FALSE;
	c1Ion				= gElementMass_x100[HYDROGEN] * 4 + gElementMass_x100[NITROGEN];

/*	Test to see if the N-terminal residue contains one amino acid, or more than one.*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(sequence[0] == gGapList[i])
		{
			oneAA = TRUE;
			break;
		}
	}
	
/*	Figure out which c1 ion to use*/
	if(oneAA)
	{
		if(sequence[1] == gGapList[Q])
		{
			c1Ion += sequence[0];
		}
		else
		{
			c1Ion = 0;	/*this is que that there was no Q*/
		}
	}
	else
	{
		testMass = sequence[0] - gGapList[Q];
		for(i = 0; i < gAminoAcidNumber; i++)
		{
			if(testMass <= gGapList[i] + gToleranceWide &&
				testMass >= gGapList[i] - gToleranceWide)
			{
				c1Ion += gGapList[i];
				break;
			}
		}
		if(c1Ion == gElementMass_x100[HYDROGEN] * 4 + gElementMass_x100[NITROGEN])
		{
			c1Ion = 0; /*que that Q is not there*/
		}
	}
	
/*	Find the c1 ion in the list of fragment masses*/
	for(i = 0; i < fragNum; i++)
	{
		if(fragMOverZ[i] >= c1Ion - gToleranceWide &&
			fragMOverZ[i] <= c1Ion + gToleranceWide)
		{
			massDiff = abs(c1Ion - fragMOverZ[i]);
			currentIonFound = ionFound[i];
			ionFound[i] = CalcIonFound(ionFound[i], massDiff);
			if(currentIonFound > ionFound[i])
			{
				ionFound[i] = currentIonFound;
			}
		}
	}
	
/*	Redefine the N-terminal 2 aa residue if it contains Q and c1*/
	if(c1Ion && !oneAA)
	{
		for(i = seqLength - 1; i > 0; i--)
		{
			sequence[i + 1] = sequence[i];
		}
		sequence[1] = gGapList[Q];
		sequence[0] = c1Ion - gElementMass_x100[HYDROGEN] * 4 - gElementMass_x100[NITROGEN];
		seqLength++;
	}
	

	return(seqLength);
}

/****************************PEFragments*********************************************
*
*	PEFragments identifies fragment ions that are due to pyridylethylated cysteines.  
*	The input is the ionFound array (0 = not identified, 1 = identified), fragNum 
*	(fragment ion number), fragmentErr (fragment ion tolerance), fragMOverZ array (of 
*	fragment ion masses), and chargeState (the precursor ion charge state).  It returns
*	nothing, but it does modify the ionFound array if ions match w/ calculated values.
*/
void PEFragments(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, 
				 INT_4 *sequence, INT_4 seqLength)
{
	INT_4 i, j, massDiff;
	INT_4 peMinusErr, pePlusErr, peFragment, peFragMinusErr, peFragPlusErr;
	REAL_4 monoCysPe;
	char test = FALSE;
	
	monoCysPe = 208.07 * gMultiplier;	/*The nominal residue mass of PE-Cys is 208.*/
	
/*	Determine if cys is present in the sequence.	*/
	for(i = 0; i < seqLength; i++)
	{
		if(sequence[i] >= (monoCysPe - gToleranceWide) &&
			sequence[i] <= (monoCysPe + gToleranceWide))	/*Check if cys is present as a single amino acid.*/
		{
			test = TRUE;
			break;
		}
	}
	if(test == FALSE)	/*Check to see if cys is present as a two amino acid gap.*/
	{
		for(i = 0; i < seqLength; i++)	/*For each "amino acid" in the array sequence.*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(sequence[i] >= (gCysPlus[j] - gToleranceWide) &&
					sequence[i] <= (gCysPlus[j] + gToleranceWide))	/*Compare each amino acid 
																		to the list of possible
																		two amino acid 
																		combinations containing 
																		cysteine.*/
				{
					test = TRUE;
					break;
				}
			}
			if(test)
			{
				break;	/*Why keep checking for cysteines if its been found already?*/
			}
		}
	}
	if(test)	/*If there is a PECys in this sequence.*/
	{
/*	Initialize variables.	*/
		peMinusErr = gParam.cysMW - gToleranceWide;
		pePlusErr = gParam.cysMW + gToleranceWide;
	
/*	Search for fragment ions.	*/
		for(i = 0; i < fragNum; i++)
		{
			if(fragMOverZ[i] <= pePlusErr)
			{
				if(fragMOverZ[i] >= peMinusErr)
				{
					massDiff = abs((long)(106.07 * gMultiplier - fragMOverZ[i]));
					ionFound[i] = CalcIonFound(ionFound[i], massDiff);
				}
			}
		}
		for(j = 1; j <= gParam.chargeState; j++)	/*Calculate the loss of PE for each charge state.*/
		{
/*	Initialize some variables.	*/
			peFragment = (gParam.peptideMW + (j * gElementMass_x100[HYDROGEN]) 
							- 105.07 * gMultiplier) / j;
			peFragMinusErr = peFragment - gToleranceWide;
			peFragPlusErr = peFragment + gToleranceWide;
		
			for(i = 0; i < fragNum; i++)	/*Find losses of PE from precursor ion.*/
			{
				if(fragMOverZ[i] >= peFragMinusErr)
				{
					if(fragMOverZ[i] <= peFragPlusErr)
					{
						massDiff = abs(peFragment - fragMOverZ[i]);
						ionFound[i] = CalcIonFound(ionFound[i], massDiff);
					}
				}
			}
		}
	}
	
	return;
}

/*****************************TotalIntensity******************************************
*
*	TotalIntensity calculates the total ion intensity found in the input array called 
*	"fragIntensity".  The region around the precursor ion is not counted in this value.  
*	It returns a INT_4 corresponding to the sum of all appropriate intensity.  Also 
*	ion intensity around the precursor ion is reassigned a value of zero.
*/

INT_4 TotalIntensity(INT_4 fragNum, INT_4 *fragMOverZ, 
						INT_4 *fragIntensity)
{
	INT_4 i;
	INT_4 totalIntensity = 0;
	INT_4 precursorMinWater, precursor, precursorMinAmmonia, highLimit;
	char charge;
	
	if(gParam.maxent3)
	{
		charge = 1;
	}
	else
	{
		charge = gParam.chargeState;
	}
	
	precursor = (gParam.peptideMW + (charge * gElementMass_x100[HYDROGEN])) / charge;
	precursorMinWater = precursor - (gWater / charge);
	precursorMinAmmonia = precursor - (gAmmonia / charge);
	highLimit = gParam.peptideMW + gElementMass_x100[HYDROGEN] - gMonoMass_x100[G] + gToleranceWide;
	
	for(i = 0; i < fragNum; i++)
	{
		if(((fragMOverZ[i] > (precursorMinWater - gToleranceWide))
			&& (fragMOverZ[i] < (precursorMinWater +  gToleranceWide))) ||
			((fragMOverZ[i] < (precursorMinAmmonia + gToleranceWide))
			&& (fragMOverZ[i] > (precursorMinAmmonia - gToleranceWide))))
		{
			fragIntensity[i] = 0;
		}
		else if(fragMOverZ[i] <= precursor + (gToleranceWide * 2)
				&& fragMOverZ[i] >= precursor - (gToleranceWide * 3))
		{
			fragIntensity[i] = 0;
		}
		else if(fragMOverZ[i] > highLimit)
		{
			fragIntensity[i] = 0;
		}
		else
		{
			totalIntensity += fragIntensity[i];
		}
	}
	
	return(totalIntensity);
}

/******************************LoadTheIonArrays********************************************
*
*	LoadTheIonArrays inputs firstMassPtr, which points to the first element in a linked 
*	list of MSData structs containing the cid data, plus two arrays defined in 
*	ScoreSequences - fragMOverZ, and fragIntensity - each of which contains MAX_ION_NUM 
*	elements.  This function initializes these two arrays.  If there are less than 
*	MAX_ION_NUM elements in the linked list of MSData structs, then these are copied 
*	directly to the two arrays.  "fragNum" is the fragment number.  If fragNum exceeds 
*	MAX_ION_NUM, then the most intense ions are loaded into the arrays.
*/

void LoadTheIonArrays(struct MSData *firstMassPtr, INT_4 *fragNum, 
						INT_4 *fragMOverZ, INT_4 *fragIntensity)
{
	struct MSData *currMSPtr, *destroyPtr;
	INT_4 i, j, destroyIndex, *mostIntMass, *mostIntInt;
	
	mostIntMass = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(mostIntMass == NULL)
	{
		printf("LoadTheIonArrays:  Out of memory.");
		exit(1);
	}
	mostIntInt = (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	if(mostIntInt == NULL)
	{
		printf("LoadTheIonArrays:  Out of memory.");
		exit(1);
	}
	
/*	
	Initialize variables.	
*/
	currMSPtr = firstMassPtr;
	*fragNum = 0;
	
/*	
	Count the number of ions in the linked list, and multiply the mass values by 100.	
*/
	while(currMSPtr != NULL)
	{
		*fragNum += 1; 
		currMSPtr = currMSPtr->next;
	}
	
/*	If the arrays are large enough, just load the linked list w/o any modifications.	*/
	if(*fragNum < MAX_ION_NUM)
	{
		currMSPtr = firstMassPtr;
		i = 0;
		while(currMSPtr != NULL)
		{
			fragMOverZ[i] = (REAL_4)currMSPtr->mOverZ + 0.5;	/*had strange effects on Mac
																where, for example, 0.5 becomes
																0.49999999 which then gets rounded
																to zero.  Only noticed when it 
																was looping thru here many times*/
			fragIntensity[i] = currMSPtr->intensity;
			i++;
			currMSPtr = currMSPtr->next;
		}
	}
	
/*	
	If there is not enough room in the arrays, then load only the most intense ions and 
	load the full array.	
*/
	else
	{
		*fragNum = MAX_ION_NUM;
		for(i = 0; i < MAX_ION_NUM; i++)		/*Find the most intense ions.*/
		{
			mostIntMass[i] = (REAL_4)firstMassPtr->mOverZ + 0.5;
			mostIntInt[i] = firstMassPtr->intensity;
			currMSPtr = firstMassPtr->next;
			while(currMSPtr != NULL)
			{
				if(currMSPtr->intensity > mostIntInt[i])
				{
					mostIntInt[i] = currMSPtr->intensity;
					mostIntMass[i] = (REAL_4)currMSPtr->mOverZ + 0.5;
					destroyPtr = currMSPtr;
				}
				currMSPtr = currMSPtr->next;
			}
			destroyPtr->intensity = 0;
		}
		for(i = 0; i < MAX_ION_NUM; i++)		/*Then sort by m/z.*/
		{
			destroyIndex = 0;
			fragMOverZ[i] = mostIntMass[0];
			fragIntensity[i] = mostIntInt[0];
			for(j = 0; j < MAX_ION_NUM; j++)
			{
				if(mostIntMass[j] < fragMOverZ[i])
				{
					fragMOverZ[i] = mostIntMass[j];
					fragIntensity[i] = mostIntInt[j];
					destroyIndex = j;
				}
			}
			mostIntMass[destroyIndex] = 10000;
		}
	}
			
	free(mostIntMass);
	free(mostIntInt);

	return;
}
		
/***************************ScoreSequences************************************
*
*	ScoreSequences is called by main in the file LutefiskMain.c and returns a linked list 
*	of structs containing information on the sequence, its rank, intensity score, and 
*	cross-correlation score.
*  
*	Input parameters are pointers to the first Sequence struct (containing the sequences 
*	to be scored), and the first MSData struct (containing the cid data), and various
*	fields within the global struct gParam - "peptideMW" (the peptide molecular weight), 
*	"fragmentErr" (the fragment ion m/z tolerance), "chargeState" (the charge state of the 
*	precursor ion), and "cysMW" (the mass of cysteine - used to figure out the type of 
*   cysteine fragmentation).  The INT_4 "rankedSeqNum" is the number of sequences that 
*   will be in the returned linked list of ranked and scored sequences (firstScorePtr), 
*   BUT IT'S STILL NOT USED.
*/
							
struct Sequence *ScoreSequences(struct Sequence *firstSequencePtr, 
						struct MSData *firstMassPtr)
{
	char cysPE, addSequence, databaseSeq;
	char argPresent = 0;
	INT_4 *sequence, *charSequence, *ionType;
	INT_4 precursor, averageBYScore = 0, highestBYScore = 0;
	INT_4 *fragMOverZ, *fragIntensity, *saveFragMOverZ;
	INT_4 i, j, fragNum, intensityTotal, seqLength = 0, storedSeqNum, countTheSeqs;
	INT_4 cleavageSites = 0, length = 0;
	REAL_4 realSeqLengthNoFudgingAtAll = 0;
	REAL_4 intScore, intOnlyScore, *ionFound, *ionFoundTemplate, *yFound, *bFound;
	REAL_4 lowMassIonConversion, lowMassCys, residueNumGuess = 0, minQuality = 0;
	REAL_8 *byError, perfectProbScore = 0;
	REAL_4 stDevErr, realSeqLength, calFactor = 1, quality = 0;
	REAL_4 probScore = 0;
	/*INT_4 m;*/	/*debug*/
	BOOLEAN aSequenceFound = FALSE; /*debugging*/
	BOOLEAN test; /*debugging*/
	INT_4 z = 0; /*debugging*/
	
	struct SequenceScore *firstScorePtr, *lowScorePtr;
	struct SequenceScore *massagedSeqListPtr = NULL, *currMassagePtr = NULL;
	struct Sequence *currSeqPtr;
	
	
/*	
*	lowMassIons contains values corresponding to amino acid immonium ions or other pieces 
*	of amino acids.
*	There are three ions per amino acid, most amino acids have a single immonium ion m/z.  
*/

	INT_4 lowMassIons[AMINO_ACID_NUMBER][3] = {
	 /* A */     440500,     0,     0, 
	 /* R */	 700657,  	870922, 1120875,  
	 /* N */	 870558,     0,     0,  
	 /* D */	 880399,     0,     0, 
	 /* C */		0, 		 0, 	0, 
	 /* E */	1020555,     0,     0,  
	 /* Q */	 840450, 1010715, 1290664,
	 /* G */	    0,       0,     0,
	 /* H */	1100718,     0,     0,  
	 /* I */	 860970, 1200483,     0,  /*I position also represent oxidized Met for qtof*/
	 /* L */	 860970,     0,     0,
	 /* K */     840814, 1011079, 1291028, 
	 /* M */	1040534,     0,     0, 
	 /* F */	1200813,     0,     0,  
	 /* P */	 700657,     0,     0,  
	 /* S */	 600449,     0,     0,
	 /* T */	 740606,     0,     0, 
	 /* W */	1590922,     0,     0, 
	 /* Y */	1360762,     0,     0,  
	 /* V */	 720813,     0,     0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0,
	 					0,0,0
	};
	
	/*check that the sequences do not exceed the array*/
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		if(currSeqPtr->peptideLength > MAX_PEPTIDE_LENGTH)
		{
			printf("LutefiskScore: peptideLength too long");
			exit(1);
		}
		currSeqPtr = currSeqPtr->next;
	}
	
/*Convert the lowMassIons to the appropriate value given the gMultiplier value.*/
	
	lowMassIonConversion = (REAL_4)gMultiplier / 10000;
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			lowMassIons[i][j] = (REAL_4)lowMassIons[i][j] * lowMassIonConversion + 0.5;
		}
	}
	
	/*Assign a value to Cys, using gParam.cysMW to account for alkylations*/
	lowMassCys = ((REAL_4)gParam.cysMW / gMultiplier);	/*convert to real mass value temporarily*/
	lowMassCys = lowMassCys - 26.9871;	/*this is loss of CO plus H*/
	lowMassCys = lowMassCys * gMultiplier;	/*convert back*/
	lowMassIons[C][0] = lowMassCys + 0.5;	/*round it off*/
	

/*
	Assign some space for the arrays.
*/	
	sequence 			= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	charSequence 		= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	fragMOverZ 			= (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	saveFragMOverZ 		= (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	fragIntensity 		= (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	ionFound 			= (float *) malloc(MAX_ION_NUM * sizeof(REAL_4));
	ionType				= (int *) malloc(MAX_ION_NUM * sizeof(INT_4));
	byError 			= (double *) malloc(MAX_ION_NUM * sizeof(REAL_8));
	ionFoundTemplate 	= (float *) malloc(MAX_ION_NUM * sizeof(REAL_4));
	bFound 				= (float *) malloc(MAX_ION_NUM * sizeof(REAL_4));
	yFound 				= (float *) malloc(MAX_ION_NUM * sizeof(REAL_4));
	
	

	
	if(gAmIHere)	/*debugging*/
	{
		aSequenceFound = CheckItOut(firstSequencePtr);
	}

/*	
	Initialize some more variables.	
*/
	
	/*For high charge states, I won't recalibrate, so use wide tolerance.*/
	if(gParam.qtofErr != 0)
	{
		if(gParam.chargeState > 3)
		{
			gParam.qtofErr = gParam.fragmentErr;
		}
	}
	gToleranceNarrow = gParam.fragmentErr * 0.95;		/*0.8This is done for the fuzzy logic.*/
	gToleranceWide = gParam.fragmentErr * 1.5;	/*1.2This is done for the fuzzy logic.*/
	storedSeqNum = 0;
	firstScorePtr = NULL;	/*This is the pointer that is returned by LutefiskScore*/
	lowScorePtr = NULL;
	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	fragNum = 0;
	gGapListDipeptideIndex = gGapListIndex;
	
	
/*	
	Check the arrays, count the sequences in the final list of
	completed sequences, put the sequence tag back into the list of sequences, determine
	if there is a C-terminal Lys or Arg for tryptic peptides, and multiply several mass
	variables by 100 so that integers can be used rather than REAL_4s.
*/
	firstSequencePtr = InitLutefiskScore(sequence, fragMOverZ, fragIntensity, ionFound, ionFoundTemplate, 
						&countTheSeqs, firstSequencePtr, yFound, bFound, firstMassPtr, charSequence,
						byError, ionType);
	
	if(gAmIHere)
	{
		aSequenceFound = CheckItOut(firstSequencePtr);
	}
	
/*
*	Figure out minimum spectrum quality based on the max sequence tag length.
*/
	residueNumGuess = gParam.peptideMW / (AV_RESIDUE_MASS * gMultiplier);
	if(residueNumGuess != 0)
	{
		minQuality = gTagLength / residueNumGuess;
	}
	
/*
*	Since the order (forward or backward) orientation of the peptide sequence is usually based
*	on the presence of characteristic ions (eg, y1=147,175 for tryptic peptides), and since
*	these ions are not necessarily very intense, and since the scores are based on accounting
*	for the greatest percentage of ion intensity, this function boosts the ion intensities
*	for these characteristic ions, if they are present.
*/
	if(gFirstTimeThru)	/*do this once; don't keep doing it for each peptide MW*/
	{
		BoostTheCTerminals(firstMassPtr);
	}
	
/*	
	Load the m/z and intensity arrays containing the cid data; count the ions, too.	
*/
	if(firstMassPtr == NULL)
	{
		printf("Most distressing!  There seems to be no CID data.");
		exit(1);
	}
	LoadTheIonArrays(firstMassPtr, &fragNum, fragMOverZ, fragIntensity);
	
/*
*	Find the average ion intensity and the standard deviation of the ion intensity,
*	and if any ions have exceed a 1.64 x stDev these intensities are adjusted
*	up or down.
*/
	
	AdjustIonIntensity(fragNum, fragIntensity);	/*since this is not affecting firstMassPtr, this can be done
												repeatedly w/o worrying about getting weird effects*/
	
	
/*	
	Initialize ionFoundTemplate elements to zero.	
*/
	for(i = 0; i < fragNum; i++)
	{
		ionFoundTemplate[i] = 0;
	}
	
/*	
	Calculate the total ion intensity. 	
*/
	intensityTotal = TotalIntensity(fragNum, fragMOverZ, fragIntensity);
									
/*	
	Figure out if the peptide is pyridylethylated.	
*/
	if((gParam.cysMW >= ((208.07 * gMultiplier) - gToleranceWide)) && 
		(gParam.cysMW <= ((208.07 * gMultiplier) + gToleranceWide)))
	{
		cysPE = 1;
	}
	else
	{
		cysPE = 0;
	}
	
/*	
	Identify losses of water or ammonia, 2 waters, 2 ammonias, or one of each from the 
	precursor ion.	
*/
	WaterLoss(ionFoundTemplate, fragNum, fragMOverZ, fragIntensity);
	
/*
	I'll first just score for b and y ions.  Those sequences that 
	seem to have alot of b and y ions delineating their sequences will be kept, whereas
	those with few b and y (or alternating b and y) will be discarded.  The first 
	sequence in the linked list will always be kept in order to
	keep things simple.
*/

	if((gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q' 
		|| gParam.fragmentPattern == 'L') && gParam.chargeState > 1)	/*Only do this for tryptic
																peptides that have precursor 
																charges greater than one.*/
	{
		if(countTheSeqs > 100)	/*don't bother weeding out ridiculous sequences if there's only a few*/
		{
			TossTheLosers(firstSequencePtr, ionFoundTemplate, fragNum, fragMOverZ, 
							fragIntensity, intensityTotal, ionFound, &countTheSeqs, 
							sequence);
		}
	
		if(gAmIHere)
		{
			aSequenceFound = CheckItOut(firstSequencePtr);
		}
							
		/*Get rid of sequences that cannot account for most of the higher m/z fragments*/
		if(countTheSeqs > 100)	/*don't bother weeding out ridiculous sequences if there's only a few*/
		{
			HighMOverZFilter(firstSequencePtr, fragMOverZ, fragIntensity, &countTheSeqs, 
							sequence, fragNum);
		}
		
		if(gAmIHere)
		{
			aSequenceFound = CheckItOut(firstSequencePtr);
		}
	}
	

/*
*	Add the database sequence(s) to the list in firstSequencePtr.
*/
	if (strlen(gParam.databaseSequences) > 0 && gCorrectMass)
	{
	    AddDatabaseSequences(firstSequencePtr);
	}
	if(gAmIHere)
	{
		aSequenceFound = CheckItOut(firstSequencePtr);
	}
	
/*
*	Yank out sequences where a two aa extension matches single aa extensions in another seq.
*	These are usually cases where a two aa extension happens to match with three amino acids.
*/

	firstSequencePtr = RemoveRedundantSequences(firstSequencePtr);
	
	if(gAmIHere)
	{
		aSequenceFound = CheckItOut(firstSequencePtr);
	}
	
/*
	For Qtof data, if the qtof error value is sufficient to distinguish between Q/K, F/M-O, and isobaric
	dipeptides, then the number of sequences are expanded to account for all of the possibilities.  Also,
	gGapList is redone to reflect the tighter qtof final score tolerance.  Since I/L are isomeric, and since
	I use L for either I/L, the 'I' position in gGapList is reserved for oxidized Met and the single letter
	code for this position is changed to 'm'.
*/
	if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0)
	{
		MakeNewgGapList();	/*Creates single aa and dipeptides and eliminates only those with exactly the same
							mass.  If any amino acids in the sequence lists are not from
							single or dipeptides (ie, three amino acids, say) then these
							masses are added to gGapList.*/
							
		/*Count the sequences, and if too many do an extensive scoring 
		using the wide tolerances of gParam.fragmentErr and remove sequences on that basis*/
		countTheSeqs = 0;
		currSeqPtr = firstSequencePtr;
		while(currSeqPtr != 0)
		{
			countTheSeqs++;
			currSeqPtr = currSeqPtr->next;
		}
		if(countTheSeqs > MAX_QTOF_SEQUENCES)
		{
			RescoreAndPrune(firstSequencePtr, ionFound, fragNum, fragMOverZ, sequence, 
		                            seqLength, argPresent, yFound, bFound, byError, 
		                            cleavageSites, lowMassIons, ionFoundTemplate,
		                            fragIntensity, intensityTotal, ionType);
		}
		if(countTheSeqs > 50 && !gCorrectMass)	/*start w/ fewer sequences if mass is wrong*/
		{
			RescoreAndPrune(firstSequencePtr, ionFound, fragNum, fragMOverZ, sequence, 
		                            seqLength, argPresent, yFound, bFound, byError, 
		                            cleavageSites, lowMassIons, ionFoundTemplate,
		                            fragIntensity, intensityTotal, ionType);
		}

		ExpandSequences(firstSequencePtr);	/*Expands the number of sequences if the qtofErr can differentiate
											between aa's and dipeptides of the same nominal mass.*/
		
		gToleranceNarrow = gParam.qtofErr * 0.95;		/*0.9This is done for the fuzzy logic.*/
		gToleranceWide = gParam.qtofErr * 3;			/*2This is done for the fuzzy logic.*/
	
	}
	
	

	/*printf("Prior to AddToGapList\n");*/
	AddToGapList(firstSequencePtr);	/*add unusual dipeptide masses to gGapList*/
	/*printf("After AddToGapList\n");*/

/*=====================================================================================
*	Here's the giant while loop, that scores each sequence in the linked list of 
*	SequenceData structs.
*/
	/*m = 0;*/	/*debug*/
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		/*m++;*/	/*debug*/
		/*if(z == 347) /*debugging*/
		/*{
			z++;
		}*/
		
/*	
*	Write the sequence to the INT_4 array "sequence", and count the number of amino acids 
*	("seqLength").  Actually, its not the sequence that is in "sequence", rather its the 
*	nominal mass of a single amino acid or a pair of amino acids times 100 (all mass values
*	are 100 x actual value in this file).
*/
		
		LoadSequence(sequence, &seqLength, currSeqPtr);
		
		/*debugging*/
		/*test = TRUE;	
		if(seqLength == 12)
		{
			for(i = 0; i < seqLength; i++)
			{
				if(gRightSequence[i] <= sequence[i] - gToleranceWide ||
					gRightSequence[i] >= sequence[i] + gToleranceWide)
				{
					test = FALSE;
				}
			}
			if(test)
			{
				i++;
			}
		}*/
		
/*	Note if this is a database sequence (used to flag such sequences in the output)*/
		
		if(currSeqPtr->gapNum == -100)
		{
			databaseSeq = TRUE;
		}
		else
		{
			databaseSeq = FALSE;
		}
		
		
/*
	Qtof data is recalibrated for each sequence using y ions greater than m/z 500.  The original fragMOverZ
	values are saved in saveFragMOverZ, and restored after the sequence has been scored.  Once the data has been
	recalibrated the tolerances are narrowed to the scorerror values from the .param file.
*/
		if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0 && gParam.chargeState <= 3)
		{
			for(i = 0; i < fragNum; i++)	/*save the original data*/
			{
				saveFragMOverZ[i] = fragMOverZ[i];
			}
			/*printf("%d Before Recalibrate\n", m);*/
			calFactor = Recalibrate(fragNum, fragMOverZ, sequence, seqLength, 
									fragIntensity);
			/*printf("%d After Recalibrate\n", m);*/
		}
		
/*	
	Initialize this variable each time around.
*/
	
		for(i = 0; i < fragNum; i++)
		{
			ionFound[i] = ionFoundTemplate[i];	/*ionFoundTemplate contains identifications 
												of ions
												that do not change from sequence to sequence.  
												Therefore, these identifications were taken 
												out of the
												giant while loop, and instead of initializing 
												"ionFound"
												to zero, its initialized to ionFoundTemplate.*/
			yFound[i] = 0;	/*yFound and bFound are used to cut the scores for sequences that
							utilze the same ions for y and b's.*/
			bFound[i] = 0;
			
			byError[i] = 100;	/*Errors in b and y ions are placed in this array, where index
								matches the fragNum indexing*/
								
			if(ionFoundTemplate[i] != 0)
			{
				ionType[i] = 1;	/*Precursor ions are of type 1*/
			}
			else
			{
				ionType[i] = 0;	/*Initialize all other positions as 0 (random) ion types*/
			}
		} 
/*	
	Identify Arg-related ions, if the precursor is singly-charged and arginine is present 
	in the sequence.  Also checks if arg is present (TRUE or FALSE).
*/

		argPresent = ArgIons(ionFound, fragNum, fragMOverZ, sequence, seqLength);

		
/* 	
	Identify a, b, and y ions.  "cleavageSites" is used in assigning the intensity based 
	score later on.  It is the number
	of times a peptide bond was cleaved or delineated by a b or y type ion.
*/
		cleavageSites = FindABYIons(ionFound, fragNum, fragMOverZ, sequence, 
		                            seqLength, argPresent, yFound, bFound, byError, ionType);
/*
	Fool with the cleavageSites value and the ionFound values for sequences that used the
	same series of ions for both y and b.
*/
		cleavageSites = AlterIonFound(ionFound, fragNum, fragMOverZ, sequence, 
		                            seqLength, yFound, bFound, cleavageSites); 
		                            
		/*if gapNum = -100, this signals that the sequence is from the database*/
		if(currSeqPtr->gapNum == -100)
		{
			cleavageSites = (currSeqPtr->peptideLength) - 1;
		}
		
/*
	Find any c1 ions when Q is at the second position from the N-terminus.
*/

		seqLength = ScoreC1(ionFound, fragNum, fragMOverZ, sequence, seqLength);
		                            
		
/*	Identify pyridylethylated cysteine fragments, but only if PE'ed cysteine is in 
	sequence.	*/

		if(cysPE)
		{
			PEFragments(ionFound, fragNum, fragMOverZ, sequence, seqLength);
		}
				
/*	Identify internal fragment ions. */

		InternalFrag(ionFound, fragMOverZ, sequence, seqLength, fragNum, ionType);

/*	Identify low mass amino acid - specific ions. */
		ScoreLowMassIons(ionFound, fragMOverZ, sequence, seqLength, lowMassIons, ionType);
	
/*	For ions that are one dalton higher than another that has been identified as b or y,
	assign this as being partially found, too.*/
		if(gParam.peakWidth < 0.7 * gMultiplier && !gParam.maxent3)
		{
			ScoreBYIsotopes(ionFound, fragMOverZ, fragNum, ionType);
		}
/*
	Calculate the actual number of amino acids in the sequence where gaps are counted as two amino acids.
*/

		realSeqLength = SequenceLengthCalc(sequence, seqLength);
		
/*
*	Calculate the number of amino acids, not accounting for Pro mis-cleavages and N-terminal dipeptides.
*/

		realSeqLengthNoFudgingAtAll = SequenceLengthCalcNoFudge(sequence, seqLength);
								
/*	Assign the intensity-based score for the sequence.*/
		if(databaseSeq)
		{
			cleavageSites = realSeqLength - 1;	/*don't penalize database-derived sequences in the score*/
		}
		intScore = IntensityScorer(fragIntensity, ionFound, cleavageSites, fragNum, 
						realSeqLength, intensityTotal);
						
/*	For qtof data, the intScore is modified so that larger correction factors to the
	calibration will attenuate the score.*/
	/*	if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0 && gParam.chargeState <= 2)
		{
			intScore = ScoreAttenuationFromCalfactor(calFactor, intScore);
		}		*/
						
/*	Determine the standard deviation of the average error.*/
		stDevErr = StandardDeviationOfTheBYErrors(byError, fragNum);
		
/*	Recalculate the intensity-based score w/o attenuation and other tricks.*/
		intOnlyScore = IntensityOnlyScorer(fragIntensity, ionFound, fragNum, intensityTotal);
/*
*	Calculate quality as mass of amino acids defined by contiguous series divided by total residue mass
*/
		quality = MassBasedQuality(sequence, seqLength, fragNum, fragMOverZ, argPresent);
		
/*
*	Determine Pavel probability score.
*/
	if(sequence[0] == 2029)
	{
		i++;
	}
		probScore = LutefiskProbScorer(sequence, seqLength, fragNum, fragMOverZ, argPresent);
		
		/*normalize the score by using the sequence length*/
		probScore = probScore / (2 * realSeqLengthNoFudgingAtAll);
		gProbScoreMax = gProbScoreMax / (2 * realSeqLengthNoFudgingAtAll);
		
		if(gProbScoreMax > perfectProbScore)
		{
			perfectProbScore = gProbScoreMax;
		}
		
/*	Qtof data is recalibrated for each sequence using y ions greater than m/z 500.  The original fragMOverZ
	values are saved in saveFragMOverZ, and restored after the sequence has been scored.  Once the data has been
	recalibrated the tolerances are narrowed to the scorerror values from the .param file.
*/
		if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0 && gParam.chargeState <= 3)
		{
			for(i = 0; i < fragNum; i++)	/*save the original data*/
			{
				fragMOverZ[i] = saveFragMOverZ[i];
			}
		}

/*	
	Check for quality of spectrum.  gCleavageSiteStringent is the longest contiguous series of b or y ions
	but it counts two aa residues the same as one aa residues.  realSeqLength is the actual lenght of the
	sequence where two aa residues are counted as two amino acids rather than a single residue. 
	gSingleAACleavageSites is the number of contiguous b or y ions that define a sequence made up entirely
	of single amino acid residues -- no gaps are included.  The realSeqLengthNoFudgingAtAll value represents
	the best guess at the number of amino acids in the sequence; in contrast, realSeqLength makes allowances
	for gaps that might include a proline, as well as an N-terminal gap -- those are counted as single residues.
*/
		/*avoid divide by zero*/
		/*if(realSeqLengthNoFudgingAtAll > 1)	
		{
				quality = gSingleAACleavageSites / (realSeqLengthNoFudgingAtAll - 1);
				length = gSingleAACleavageSites;
		}
		else
		{
			quality = 0;
			length = 0;
		}		This is the length based quality, which is not valid unless unsequenced regions are dipeptides*/
		
		
		
/*
*	Adjust quality value higher if its zero and a sequence tag was found
*/	

		if(quality == 0)
		{
			if(minQuality > 0)
			{
				quality = minQuality;	/*minQuality is derived from the fraction of sequence covered
										by a sequence tag*/
			}
		}
		
/*	Store the sequence, intensity score, actual peptide length, and quality.*/

		addSequence = IsThisADuplicate(firstScorePtr, sequence, intOnlyScore, 
									intScore, seqLength);
		
		if(addSequence || currSeqPtr->gapNum == -100)	/*-100 flag for database seq*/
		{
			if(storedSeqNum <= MAX_X_CORR_NUM) 
			{
	
				firstScorePtr = AddToSeqScoreList(firstScorePtr, LoadSeqScoreStruct(intScore,
								intOnlyScore, sequence, charSequence, seqLength, stDevErr, 
								cleavageSites, calFactor, databaseSeq, intOnlyScore, quality, length,
								probScore, 0.0));
				storedSeqNum++;
	
			}
			else
			{
				if(lowScorePtr == NULL)  /*Find the lowest intensity-based score out of all 
										stored sequences.*/
				{
					lowScorePtr = FindLowestScore(firstScorePtr);
				}
				if(intScore > lowScorePtr->intensityScore)
				{
					lowScorePtr->intensityScore = intScore;
					for(i = 0; i < seqLength; i++)
					{
						lowScorePtr->peptide[i] = sequence[i];
						lowScorePtr->peptideSequence[i] = charSequence[i];
					}
					lowScorePtr->peptide[seqLength] = 0;
					lowScorePtr->intensityOnlyScore = intOnlyScore;
					lowScorePtr->stDevErr = stDevErr;
					lowScorePtr->crossDressingScore = intScore;
					lowScorePtr->calFactor = calFactor;
					lowScorePtr->quality = quality;
					lowScorePtr->length = length;
					lowScorePtr->probScore = probScore;
					lowScorePtr->cleavageSites = cleavageSites;
					lowScorePtr->databaseSeq = databaseSeq;
					lowScorePtr->rank = 0;
					lowScorePtr = NULL;	/*If this is not NULL, then that means it found the lowScorePtr
										earlier, but the sequence that was previously under consideration
										had a lower score.  This means keeps the program from searching
										for the same lowScorePtr, unless its been NULLed.*/
				}
			}
		}
		
		/*debugging*/
		if(gAmIHere)
		{
			aSequenceFound = CheckItOutSequenceScore(firstScorePtr);
			if(!aSequenceFound)
			{
				z++;	/*stop in debugger*/
			}
			else
			{
				z++;
			}
		}/*end debugging*/
				
		currSeqPtr = currSeqPtr->next;	/*Point to the next struct in the linked list to 
										continue the giant while loop.*/
	}	/*End of the giant while loop.*/
	
	if(gAmIHere)
	{
		aSequenceFound = CheckItOutSequenceScore(firstScorePtr);
	}

/*
*	Convert relevant mass values back to the real masses (not the * gMultiplier values).  Convert
*	the mass-based peptide sequence to a character-based sequence and place in the peptideSequence 
*	field of firstScorePtr.
*/

	RevertBackToReals(firstMassPtr, firstScorePtr);

	if(gAmIHere)
	{
		aSequenceFound = CheckItOutSequenceScore(firstScorePtr);
	}

/*	
*	Here's where a rank is assigned based on the intensity-based score.
*/

	SeqIntensityRanker(firstScorePtr);
	
	if(gAmIHere)
	{
		aSequenceFound = CheckItOutSequenceScore(firstScorePtr);
	}
	
/*	
*	Here's where the cross-correlation scoring would be done - after the giant while loop.
*	This is where the top rankedSeqNum number of scores are found, and the rest are discarded.  
*	Cross-correlation scores are normalized to the auto-correlation of the actual spectrum.  The
*	background associated with tau = 0 is determined by adding the symetrical differences around tau = 0.
*	This takes advantage of the fact that a real match is going to be symetrical, and bogus matches are not.
*/
	if (gParam.fMonitor && gCorrectMass)
	{
		printf("Cross-dressing.\n");
		fflush(stdout);
	}
	
	DoCrossCorrelationScoring(firstScorePtr, firstMassPtr);
	
/*
*	Figure out a theoretically perfect probScore for comparison with actual probScores.
*/
	
	/*perfectProbScore = CalcPerfectProbScore(fragNum, fragMOverZ);*/
	
	if(gAmIHere)
	{
		aSequenceFound = CheckItOutSequenceScore(firstScorePtr);
	}
	
	if(gParam.fMonitor && gCorrectMass)
	{
		PrintToConsole(firstScorePtr);
				
/*		PrintScoreDetailsToXLFile(firstScorePtr, perfectProbScore);*/	/*debugging output*/
	}
	
	
/*	Massage the scores to come up with a short list of most likely sequences. */

/*	massagedSeqListPtr = MassageScores(firstScorePtr);*/
	
	massagedSeqListPtr = DetermineBestCandidates(firstScorePtr);
	
/*
	Find max quality and length values from the final list in the output.
*/
	currMassagePtr = massagedSeqListPtr;
	quality = 0;
	length = 0;
	while(currMassagePtr != NULL)
	{
		if(currMassagePtr->quality > quality)
		{
			quality = currMassagePtr->quality;
			length = currMassagePtr->length;
		}
		currMassagePtr = currMassagePtr->next;
	}

/*  Stop the clock */
	if(gCorrectMass)
	{
		gParam.searchTime = (clock() - gParam.startTicks)/ CLOCKS_PER_SEC;
	}

/*	Output is printed to the console and to a file.*/
	if(gCorrectMass)
	{
		PrintToConsoleAndFile(massagedSeqListPtr, quality, length, perfectProbScore);
		
	}
	
/*	Free some linked lists*/
	FreeAllSequenceScore(firstScorePtr);
	FreeAllSequenceScore(massagedSeqListPtr);
	
	RevertTheRevertBackToReals(firstMassPtr);
	
/*	Free the arrays, before I forget.*/
	free(sequence);
	free(fragMOverZ);
	free(fragIntensity);
	free(ionFound);	
	free(ionFoundTemplate);	
	free(yFound);
	free(bFound);
	free(byError);
	free(charSequence);
	free(saveFragMOverZ);

	return(firstSequencePtr);		/*Return a pointer to the massaged list of sequences and scores.*/
}
					
/***************************MassBasedQuality**********************************
*
*	Calculates quality as the mass of amino acids defined by a contiguous ion
*	series, divided by total mass of all residues.
*/
REAL_4	MassBasedQuality(INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *fragMOverZ, char argPresent)
{
	INT_4 i, j, k, maxCharge;
	REAL_4 bQualMass, yQualMass, quality, bSinglyCharged, ySinglyCharged, bIon, yIon, totalResidueMass;
	REAL_4 maxBQualMass, maxYQualMass, lowMassLimit;
	REAL_4 bMinWater, bMinAmmonia, yMinWater, yMinAmmonia;
	BOOLEAN singleAA, bIonFound, yIonFound, validTerminus;
	
	totalResidueMass = 0;
	maxBQualMass = 0;
	maxYQualMass = 0;
	lowMassLimit = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;
	lowMassLimit = lowMassLimit * 0.33333;
	if(gParam.chargeState > 1)
	{
		maxCharge = gParam.chargeState - 1;
	}
	else
	{
		maxCharge = 1;
	}
	
	/*Calculate quality based on b ions*/
	bQualMass = 0;
	bSinglyCharged = gParam.modifiedNTerm;
	
	for(i = 0; i < seqLength; i++)
	{
		singleAA = FALSE;
		validTerminus = FALSE;
		bSinglyCharged = bSinglyCharged + sequence[i];
		totalResidueMass += sequence[i];
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] <= gMonoMass_x100[j] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] - gToleranceWide)
			{
				singleAA = TRUE;
				validTerminus = TRUE;
				break;
			}
		}
		
		bIonFound = FALSE;	/*set it to FALSE for each "residue"*/
		
		if(i == 0 && !singleAA)	/*Count the two aa at the N-terminus*/
		{
			for(j = 0; j < gGapListDipeptideIndex; j++)
			{
				if(sequence[0] == gGapList[j])	/*make sure its really two aa, and not something more*/
				{
					bIonFound = TRUE;
					validTerminus = TRUE;
					break;
				}
			}
		}
		
		
		if(singleAA)
		{
			for(j = 1; j <= maxCharge; j++)
			{
				bIon = (bSinglyCharged + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				bMinWater = bSinglyCharged - (gElementMass_x100[HYDROGEN] * 2 + gElementMass_x100[OXYGEN]);
				bMinAmmonia = bSinglyCharged - (gElementMass_x100[HYDROGEN] * 3 + gElementMass_x100[NITROGEN]);
				bMinWater = (bMinWater + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				bMinAmmonia = (bMinAmmonia + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				for(k = 0; k < fragNum; k++)
				{
					if(bIon <= fragMOverZ[k] + gToleranceWide && bIon >= fragMOverZ[k] - gToleranceWide)
					{
						bIonFound = TRUE;
						break;
					}
					if(argPresent)
					{
						if(bMinWater <= fragMOverZ[k] + gToleranceWide 
							&& bMinWater >= fragMOverZ[k] - gToleranceWide)
						{
							bIonFound = TRUE;
							break;
						}
						if(bMinAmmonia <= fragMOverZ[k] + gToleranceWide 
							&& bMinAmmonia >= fragMOverZ[k] - gToleranceWide)
						{
							bIonFound = TRUE;
							break;
						}
					}
				}
				if(bIonFound)
				{
					break;
				}
			}
			/*what if its ion trap data and the 1/3 rule predicts that the ions are too low of mass to be seen?*/
			if(!bIonFound && gParam.fragmentPattern == 'L' && bSinglyCharged < lowMassLimit)
			{
				for(j = 1; j <= maxCharge; j++)
				{
					yIon = gParam.peptideMW - bSinglyCharged + 2 * gElementMass_x100[HYDROGEN];
					yMinWater = yIon - (gElementMass_x100[OXYGEN] + 2 * gElementMass_x100[HYDROGEN]);
					yMinAmmonia = yIon - (gElementMass_x100[NITROGEN] + 3 * gElementMass_x100[HYDROGEN]);
					
					yIon = (yIon + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					yMinWater = (yMinWater + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					yMinAmmonia = (yMinAmmonia + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					
					for(k = 0; k < fragNum; k++)
					{
						if(yIon <= fragMOverZ[k] + gToleranceWide && yIon >= fragMOverZ[k] - gToleranceWide)
						{
							bIonFound = TRUE;
							break;
						}
						if(argPresent)
						{
							if(yMinWater <= fragMOverZ[k] + gToleranceWide && 
								yMinWater >= fragMOverZ[k] - gToleranceWide)
							{
								bIonFound = TRUE;
								break;
							}
							if(yMinAmmonia <= fragMOverZ[k] + gToleranceWide && 
								yMinAmmonia >= fragMOverZ[k] - gToleranceWide)
							{
								bIonFound = TRUE;
								break;
							}
						}
					}
				}
			}
		}
		if(bIonFound || (validTerminus && i == seqLength-1)) /*i==seqLenght-1 is the b ion of MH+-18*/
		{
			bQualMass = bQualMass + sequence[i];
		}
		else
		{
			if(bQualMass > maxBQualMass)
			{
				maxBQualMass = bQualMass;
			}
			bQualMass = 0;
		}
	}

	
	/*Calculate quality based on y ions*/
	yQualMass = 0;
	ySinglyCharged = gParam.modifiedCTerm + 2 * gElementMass_x100[HYDROGEN];
	
	for(i = seqLength - 1; i >= 0; i--)
	{
		singleAA = FALSE;
		validTerminus = FALSE;
		ySinglyCharged = ySinglyCharged + sequence[i];
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] <= gMonoMass_x100[j] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] - gToleranceWide)
			{
				singleAA = TRUE;
				validTerminus = TRUE;
				break;
			}
		}
		
		yIonFound = FALSE;	/*set it to FALSE for each "residue"*/
		
		
		if(i == 0 && !singleAA)	/*Count the two aa at the N-terminus*/
		{
			for(j = 0; j < gGapListDipeptideIndex; j++)
			{
				if(sequence[0] == gGapList[j])	/*make sure its really two aa, and not something more*/
				{
					yIonFound = TRUE;
					validTerminus = TRUE;
					break;
				}
			}
		}
		
		if(singleAA)
		{
			for(j = 1; j <= maxCharge; j++)
			{
				yIon = (ySinglyCharged + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				yMinWater = ySinglyCharged - (gElementMass_x100[HYDROGEN] * 2 + gElementMass_x100[OXYGEN]);
				yMinAmmonia = ySinglyCharged - (gElementMass_x100[HYDROGEN] * 3 + gElementMass_x100[NITROGEN]);
				yMinWater = (yMinWater + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				yMinAmmonia = (yMinAmmonia + (j-1) * gElementMass_x100[HYDROGEN]) / j;
				for(k = 0; k < fragNum; k++)
				{
					if(yIon <= fragMOverZ[k] + gToleranceWide && yIon >= fragMOverZ[k] - gToleranceWide)
					{
						yIonFound = TRUE;
						break;
					}
					if(argPresent)
					{
						if(yMinWater <= fragMOverZ[k] + gToleranceWide 
							&& yMinWater >= fragMOverZ[k] - gToleranceWide)
						{
							yIonFound = TRUE;
							break;
						}
						if(yMinAmmonia <= fragMOverZ[k] + gToleranceWide 
							&& yMinAmmonia >= fragMOverZ[k] - gToleranceWide)
						{
							yIonFound = TRUE;
							break;
						}
					}
				}
				if(yIonFound)
				{
					break;
				}
			}
			/*what if its ion trap data and the 1/3 rule predicts that the ions are too low of mass to be seen?*/
			if(!yIonFound && gParam.fragmentPattern == 'L' && ySinglyCharged < lowMassLimit)
			{
				for(j = 1; j <= maxCharge; j++)
				{
					bIon = gParam.peptideMW - ySinglyCharged + 2 * gElementMass_x100[HYDROGEN];
					bMinWater = bIon - (gElementMass_x100[OXYGEN] + 2 * gElementMass_x100[HYDROGEN]);
					bMinAmmonia = bIon - (gElementMass_x100[NITROGEN] + 3 * gElementMass_x100[HYDROGEN]);
					
					bIon = (bIon + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					bMinWater = (bMinWater + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					bMinAmmonia = (bMinAmmonia + (j-1) * gElementMass_x100[HYDROGEN]) / j;
					
					for(k = 0; k < fragNum; k++)
					{
						if(bIon <= fragMOverZ[k] + gToleranceWide && bIon >= fragMOverZ[k] - gToleranceWide)
						{
							yIonFound = TRUE;
							break;
						}
						if(argPresent)
						{
							if(bMinWater <= fragMOverZ[k] + gToleranceWide && 
								bMinWater >= fragMOverZ[k] - gToleranceWide)
							{
								yIonFound = TRUE;
								break;
							}
							if(bMinAmmonia <= fragMOverZ[k] + gToleranceWide && 
								bMinAmmonia >= fragMOverZ[k] - gToleranceWide)
							{
								yIonFound = TRUE;
								break;
							}
						}
					}
				}
			}
		}
		if(yIonFound || (validTerminus && i==0))	/*i==0 means that I don't have to look for MH+*/
		{
			yQualMass = yQualMass + sequence[i];
		}
		else
		{
			if(yQualMass > maxYQualMass)
			{
				maxYQualMass = yQualMass;
			}
			yQualMass = 0;
		}
	}
	
	if(yQualMass > maxYQualMass)
	{
		maxYQualMass = yQualMass;
	}
	if(bQualMass > maxBQualMass)
	{
		maxBQualMass = bQualMass;
	}
	
	if(maxYQualMass > maxBQualMass)
	{
		quality = maxYQualMass / totalResidueMass;
	}
	else
	{
		quality = maxBQualMass / totalResidueMass;
	}

	return(quality);
}

					
/***************************RescoreAndPrune***********************************
*
*
*/
void	RescoreAndPrune(struct Sequence *firstSequencePtr, REAL_4 *ionFound, INT_4 fragNum, 
									INT_4 *fragMOverZ, INT_4 *sequence, INT_4 seqLength, 
									char argPresent, REAL_4 *yFound, REAL_4 *bFound, 
									REAL_8 *byError, INT_4 cleavageSites, 
									INT_4 lowMassIons[][3], REAL_4 *ionFoundTemplate, 
									INT_4 *fragIntensity, INT_4 intensityTotal, INT_4 *ionType)
		                         
{
struct Sequence *currSeqPtr = NULL;
struct Sequence *previousPtr = NULL;
INT_4 countTheSeqs = 0;
INT_4 i, realSeqLength;
REAL_4 maxScoreFraction = 0.45;
REAL_4 intScore, maxScore;

currSeqPtr = firstSequencePtr;
while(currSeqPtr != NULL)
{

	LoadSequence(sequence, &seqLength, currSeqPtr);
	for(i = 0; i < fragNum; i++)
	{
		ionFound[i] = ionFoundTemplate[i];
		yFound[i] = 0;	
		bFound[i] = 0;
		byError[i] = 100;	
	} 
	
	cleavageSites = FindABYIons(ionFound, fragNum, fragMOverZ, sequence, 
		                         seqLength, argPresent, yFound, bFound, byError, ionType);
		                            
	cleavageSites = AlterIonFound(ionFound, fragNum, fragMOverZ, sequence, 
		                           seqLength, yFound, bFound, cleavageSites); 
	
	/*if gapNum = -100, this signals that the sequence is from the database*/
	if(currSeqPtr->gapNum == -100)
	{
		cleavageSites = (currSeqPtr->peptideLength) - 1;
	}
	
	InternalFrag(ionFound, fragMOverZ, sequence, seqLength, fragNum, ionType);

	ScoreLowMassIons(ionFound, fragMOverZ, sequence, seqLength, lowMassIons, ionType);
	
	if(gParam.peakWidth < 0.7 * gMultiplier && !gParam.maxent3)
	{
		ScoreBYIsotopes(ionFound, fragMOverZ, fragNum, ionType);
	}

	realSeqLength = SequenceLengthCalc(sequence, seqLength);
								
	intScore = IntensityScorer(fragIntensity, ionFound, cleavageSites, fragNum, 
						realSeqLength, intensityTotal);
						
	currSeqPtr->score = intScore * 1000 + 0.5;
	
	currSeqPtr = currSeqPtr->next;
}

/*	Find max score*/
	currSeqPtr = firstSequencePtr;
	maxScore = currSeqPtr->score;
	while(currSeqPtr != NULL)
	{
		if(currSeqPtr->score > maxScore)
		{
			maxScore = currSeqPtr->score;
		}
		currSeqPtr = currSeqPtr->next;
	}
	
/*	Count the sequences*/
	countTheSeqs = 0;
	currSeqPtr = firstSequencePtr;
	while(currSeqPtr != NULL)
	{
		countTheSeqs++;
		currSeqPtr = currSeqPtr->next;
	}
	
/*	Remove sequences with less than 0.7 x maxScore*/
	while(countTheSeqs > MAX_QTOF_SEQUENCES)
	{
		maxScoreFraction = maxScoreFraction + 0.02;
		currSeqPtr = firstSequencePtr->next;
		previousPtr = firstSequencePtr;
		while(currSeqPtr != NULL)
		{
			if(currSeqPtr->score < maxScore * maxScoreFraction)
			{
				previousPtr->next = currSeqPtr->next;
				free(currSeqPtr);
				currSeqPtr = previousPtr->next;
			}
			else
			{
				currSeqPtr = currSeqPtr->next;
				previousPtr = previousPtr->next;
			}
		
		}
		/*	Count the sequences*/
		countTheSeqs = 0;
		currSeqPtr = firstSequencePtr;
		while(currSeqPtr != NULL)
		{
			countTheSeqs++;
			currSeqPtr = currSeqPtr->next;
		}
	}
	
	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Scoring %4ld sequences found for qtof scoring.\n", countTheSeqs);
		printf("These had intensity scores in excess of %.3f.\n", maxScoreFraction);
	}
		
	return;
}







/********************************LutefiskProbScorer*****************************************************
*
*	Assign probability scores to sequences.
*
*/
REAL_4		LutefiskProbScorer(INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *fragMOverZ, char argPresent)
{
	INT_4 	i;
	REAL_4	*randomProb;
	REAL_8	probScore = 0;

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
	/*debug*/
	if(sequence[0] == 7104)
	{
		i++;
	}
			
/*	Calculate random probability for each ion*/

	CalcRandomProb(randomProb, fragMOverZ, fragNum);
		
/*	Score the sequences*/
	
	/*Get initial probability based on terminal group (Lys and Arg are good; others are not)*/
	probScore = InitProbScore(sequence, seqLength);
	
	/*Initialize the maximum probability score possible for this sequence (used to normalize later)*/
	gProbScoreMax = probScore;
	
	/*Find the b ions*/
	probScore = FindBIons(fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength, argPresent);
	
	/*Find the y ions*/
	probScore = FindYIons(fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength, argPresent);
	
	/*Find the internal fragment ions*/
	probScore = FindInternalIons(fragMOverZ, fragNum, probScore, randomProb, sequence, seqLength);
	
	/*Find the immonium ions*/
	probScore = FindImmoniumIons(fragMOverZ, fragNum, probScore, randomProb,
									sequence, seqLength);
		
	/*Change probability score to log base 10 scale*/
	if(probScore > 1)
	{
		probScore = log10(probScore);
	}
	else	/*keep things positive by only logging things over a value of 1*/
	{
		probScore = 0.0001;
	}
	
	/*gProbScoreMax is based on y and b ion scoring, and assumes all reasonable values were found*/
	if(gProbScoreMax > 1)
	{
		gProbScoreMax = log10(gProbScoreMax);
	}
	else	/*keep things positive by only logging things over a value of 1*/
	{
		gProbScoreMax = .0001;
	}
	
/*	probScore = probScore / gProbScoreMax;*/
		
/*	Free array*/
	free(randomProb);
	
	return(probScore);
}


/***********************************FindImmoniumIons******************************************
*
*	Finds and scores the amino acid immonium ions.
*/
REAL_8	FindImmoniumIons(INT_4 *mass, INT_4 ionCount, 
							REAL_8 probScore, REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength)
{	
	REAL_4 lowMassIons[AMINO_ACID_NUMBER][3] = {
	 /* A */     0,			0,			0, 
	 /* R */	 0,			0,			0,  
	 /* N */	 0,			0,			0,  
	 /* D */	 0,			0,			0, 
	 /* C */	 0,			0,			0, 
	 /* E */	 0,			0,			0,  
	 /* Q */	 84.0450, 	101.0715,	129.0664,
	 /* G */	 0,			0,			0,
	 /* H */	110.0718,	0,			0,  
	 /* I */	 86.0970,	0,			0,  
	 /* L */	 86.0970,	0,			0,
	 /* K */     84.0814, 	101.1079, 	129.1028, 
	 /* M */	104.0534,	0,			0, 
	 /* F */	120.0813,	0,			0,  
	 /* P */	 70.0657,	0,			0,  
	 /* S */	 0,			0,			0,
	 /* T */	 0,			0,			0, 
	 /* W */	159.0922,	0,			0, 
	 /* Y */	136.0762,	0,			0,  
	 /* V */	 72.0813,	0,			0,
	 			 0,			0,			0,
	 			 0,			0,			0,
	 			 0,			0,			0,
	 			 0,			0,			0,
	 			 0,			0,			0
	};
	INT_4 i, j, k, immoniumIndex, aaPresent[AMINO_ACID_NUMBER];
	REAL_4 individualProb;
	REAL_4 massDiff, errProb, testErrProb;
	BOOLEAN areThereAnyLowMassIons = FALSE;
	BOOLEAN immoniumFound = FALSE;

	/*multiply low mass values by gMultiplier*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(lowMassIons[i][j] != 0)
			{
				lowMassIons[i][j] *= gMultiplier;
			}
		}
	}
	
	/*Check to see if there are any immonium ions at all*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(lowMassIons[i][j] != 0)
			{
				for(k = 0; k < ionCount; k++)
				{
					if(mass[k] > lowMassIons[W][0] + gToleranceWide)
					{
						break;
					}
					if(mass[k] >= lowMassIons[i][j] - gToleranceWide &&
						mass[k] <= lowMassIons[i][j] + gToleranceWide)
					{
						areThereAnyLowMassIons = TRUE;
						break;
					}
				}
			}
		}
	}
	
	/*Figure out which amino acids are present*/
	for(i = 0; i < AMINO_ACID_NUMBER; i++)
	{
		aaPresent[i] = 0;	
	}
	
	for(i = 0; i < seqLength; i++)
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] <= gMonoMass_x100[j] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] - gToleranceWide)
			{
				aaPresent[j] = 1;
				break;
			}
		}
	}
	
	/*Proceed if there are any immonium ions at all*/
	if(areThereAnyLowMassIons)
	{
		for(i = 0; i < gAminoAcidNumber; i++)
		{
			if(aaPresent[i] && lowMassIons[i][0] != 0)
			{
				errProb = 0;	/*use best match when several ions exist for a given amino acid*/
				for(j = 0; j < 3; j++)
				{
					if(lowMassIons[i][j] > 0)
					{
						immoniumIndex = 0;
						immoniumFound = FALSE;
						for(k = 0; k < ionCount; k++)
						{
							if(mass[k] > lowMassIons[W][0] + gToleranceWide)
							{
								break;
							}
							if(mass[k] < lowMassIons[i][j] + gToleranceWide &&
								mass[k] > lowMassIons[i][j] - gToleranceWide)
							{
								immoniumIndex = k;
								immoniumFound = TRUE;
								massDiff = fabs(lowMassIons[i][j] - mass[k]);
								testErrProb = CalcIonFound(0, massDiff);
								if(testErrProb > errProb)
								{
									errProb = testErrProb;
								}
							}
						}
					}
				}
						
				if(immoniumFound)	/*something was found*/
				{
					individualProb = immoniumProb / randomProb[immoniumIndex];
					individualProb *= errProb;
					if(individualProb > 1)
					{
						probScore *= individualProb;
					}
				}
				else	/*immonium ions not found, so penalize*/
				{
					individualProb = (1 - immoniumProb) / (1 - randomProb[immoniumIndex]);
					if(individualProb < 1 && individualProb > 0)
					{
						probScore *= individualProb;
					}
				}
			}
		}
	}

	return(probScore);
}
/***********************************FindInternalIons******************************************
*
*	Finds and scores the internal fragment ions.
*/
REAL_8	FindInternalIons(INT_4 *mass, INT_4 ionCount, 
						REAL_8 probScore, REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength)
{			
	INT_4 i, j, k, residueCount;
	REAL_4 testMass, individualProb;
	REAL_4 massDiff, errProb;
	REAL_4 precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;
	BOOLEAN nTermPro, intFragTest;
	
	if(seqLength < 4)
	{
		return(probScore);	/*need at least four residues for an internal fragment*/
	}
	
	for(i = 1; i < seqLength - 2; i++)
	{
		testMass = sequence[i] + gElementMass_x100[HYDROGEN];
		residueCount = 1;
		if(sequence[i] > gMonoMass_x100[P] - gToleranceWide &&
			sequence[i] < gMonoMass_x100[P] + gToleranceWide)
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
			if(testMass < precursor - gToleranceWide && residueCount < 5
					&& testMass > mass[0])	/*dont bother w/ high mass internal frags*/
			{
				intFragTest = FALSE;
				for(k = 0; k < ionCount; k++)
				{	
					if(mass[k] > testMass + gToleranceWide)
					{
						break;	/*I need to save the k value at the point where this occurs*/
					}
					if(testMass < mass[k] + gToleranceWide &&
						testMass > mass[k] - gToleranceWide)
					{
						intFragTest = TRUE;
						massDiff = fabs(testMass - mass[k]);
						errProb = CalcIonFound(0, massDiff);
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
						individualProb *= errProb;
						if(individualProb > 1)
						{
							probScore *= individualProb;
						}
					}
					else
					{
						individualProb = internalProb / randomProb[k];
						individualProb *= errProb;
						if(individualProb > 1)
						{
							probScore *= individualProb;
						}
					}
				}
				else	/*didn't find any, so penalize*/
				{
					if(nTermPro)
					{
						individualProb = (1 - internalProProb) / (1 - randomProb[k]);
						if(individualProb < 1 && individualProb > 0)
						{
							probScore *= individualProb;
						}
					}
					else
					{
						individualProb = (1 - internalProb) / (1 - randomProb[k]);
						if(individualProb < 1 && individualProb > 0)
						{
							probScore *= individualProb;
						}
					}
				}
			}
		}
	}
	
	return(probScore);
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
	
	
	
	residueMass = sequence[seqLength - 1];
	
	/*Initialize for tryptic peptides*/
	if(gParam.proteolysis == 'T')
	{
		if(residueMass < gMonoMass_x100[R] + gToleranceWide &&
			residueMass > gMonoMass_x100[R] - gToleranceWide)
		{
			score = 0.95;
		}
		else if(residueMass < gMonoMass_x100[K] + gToleranceWide &&
			residueMass > gMonoMass_x100[K] - gToleranceWide)
		{
			score = 0.95;
		}
		else
		{
			for(i = 0; i < gAminoAcidNumber; i++)
			{
				testMass = residueMass - gMonoMass_x100[i];
				if(testMass < gMonoMass_x100[R] + gToleranceWide &&
					testMass > gMonoMass_x100[R] - gToleranceWide)
				{
					score = 0.95;
					break;
				}
				else if(testMass < gMonoMass_x100[K] + gToleranceWide &&
					testMass > gMonoMass_x100[K] - gToleranceWide)
				{
					score = 0.95;
					break;
				}
			}
		}
	}
	else if(gParam.proteolysis == 'K')
	{
		if(residueMass < gMonoMass_x100[K] + gToleranceWide &&
			residueMass > gMonoMass_x100[K] - gToleranceWide)
		{
			score = 0.95;
		}
		else
		{
			for(i = 0; i < gAminoAcidNumber; i++)
			{
				testMass = residueMass - gMonoMass_x100[i];
				if(testMass < gMonoMass_x100[K] + gToleranceWide &&
					testMass > gMonoMass_x100[K] - gToleranceWide)
				{
					score = 0.95;
					break;
				}
			}
		}
	}
	else if(gParam.proteolysis == 'E')
	{
		if(residueMass < gMonoMass_x100[E] + gToleranceWide &&
			residueMass > gMonoMass_x100[E] - gToleranceWide)
		{
			score = 0.95;
		}
		else if(residueMass < gMonoMass_x100[D] + gToleranceWide &&
			residueMass > gMonoMass_x100[D] - gToleranceWide)
		{
			score = 0.95;
		}
		else
		{
			for(i = 0; i < gAminoAcidNumber; i++)
			{
				testMass = residueMass - gMonoMass_x100[i];
				if(testMass < gMonoMass_x100[E] + gToleranceWide &&
					testMass > gMonoMass_x100[E] - gToleranceWide)
				{
					score = 0.95;
					break;
				}
				else if(testMass < gMonoMass_x100[D] + gToleranceWide &&
					testMass > gMonoMass_x100[D] - gToleranceWide)
				{
					score = 0.95;
					break;
				}
			}
		}
	}

	return(score);
}

/************************************CalcRandomProb*********************************************
*
*	For each ion, a 400 u window is identified (usually +/- 200 u surrounding it) and the number
*	of ions is counted within the window.  That counted number is divided by the number of possible
*	ions that could fit in that 400 u window, which depends on the instrument resolution.
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
		if(mass[i] < lowMass + 200 * gMultiplier)	/*bottom 400 u window before it moves*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] < lowMass + 400 * gMultiplier)
				{
					windowCount++;
				}
			}
		}
		else if(mass[i] > highMass - 200 * gMultiplier)	/*top 400 u window that stops moving*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] > highMass - 400 * gMultiplier)
				{
					windowCount++;
				}
			}
		}
		else	/*this is the moving window*/
		{
			for(j = 0; j < ionCount; j++)
			{
				if(mass[j] > mass[i] - 200 * gMultiplier && mass[j] < mass[i] + 200 * gMultiplier)
				{
					windowCount++;
				}
			}
		}
		/*calculate the randomness of this ion*/
		randomProb[i] = (REAL_4) windowCount / 400;	/*assuming low resolution of one peak per amu is possible*/
		
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


/**************************************FindBIons*********************************************
*
*	Find the b ions for the sequence and change the ionFound to 1.  Return a value that corresponds
*	to the number of consecutive b ions.
*/

REAL_8 FindBIons(INT_4 *mass, INT_4 ionCount, REAL_8 probScore,
					REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength, char argPresent)
{	
	INT_4	i, j, k, bIonIndex, posResidues, proGapLimit, proGapCount, oxMetCount, maxCharge;
	INT_4	addBIons;
	REAL_4	water, ammonia, bIonTemplate, bIonMin17Template, bIonMin18Template, precursor;
	REAL_4	bIon, bIonMin17, bIonMin18, individualProb, aIon, aIonTemplate, carbonMonoxide;
	REAL_4	neutralLossProb, highMassLimit, lowMassLimit, bMin64;
	REAL_4	bErrProb, b17ErrProb, b18ErrProb, aErrProb, massDiff, b64ErrProb;
	REAL_4	lossOf64, bIonMin64Template;
	BOOLEAN	bIonTest, bMin18Test, bMin17Test, aIonTest, isItAGap, bMin64Test;
	BOOLEAN nTerminalQ, nTerminalE, TwoAAGap;
	
	
/*	Initialize*/
	water = gElementMass_x100[OXYGEN] + gElementMass_x100[HYDROGEN] * 2;
	ammonia = gElementMass_x100[NITROGEN] + gElementMass_x100[HYDROGEN] * 3;
	carbonMonoxide = gElementMass_x100[CARBON] + gElementMass_x100[OXYGEN];
	lossOf64 = gElementMass_x100[CARBON] + gElementMass_x100[SULFUR] + 
				gElementMass_x100[OXYGEN] + 4 * gElementMass_x100[HYDROGEN];
	bIonTemplate = gParam.modifiedNTerm;
	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;
	posResidues = 1;
	oxMetCount = 0;
	if(gParam.peptideMW < 1000 * gMultiplier)
	{
		proGapLimit = 1;	/*gaps w/ Pro are not counted as gaps, unless there are more than 1 such Pro gap*/
	}
	else
	{
		proGapLimit = 2;	/*the limit is higher for higher molecular weight peptides*/
	}
	proGapCount = 0;
	if(gParam.chargeState == 1)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = 2;
	}
	
	if(gParam.fragmentPattern == 'L')	/*ion masses outside of this range are not penalized if not found*/
	{
		lowMassLimit = precursor * 0.333;	/*so-called 1/3 rule*/
		highMassLimit = 2000 * gMultiplier;	/*mass limit for Deca*/
	}
	else
	{
		lowMassLimit = 146 * gMultiplier;	/*y1 for Lys*/
		highMassLimit = 2 * precursor;	/*often the very high mass ions are missing*/
	}
	
/*Determine if the N-terminal residue is Q or E.  If so, then b-17 and b-18 are counted even if no b*/
	if(sequence[0] >= gMonoMass_x100[Q] - gToleranceWide 
		&& sequence[0] <= gMonoMass_x100[Q] + gToleranceWide)
	{
		nTerminalQ = TRUE;
	}
	else
	{
		nTerminalQ = FALSE;
	}
	if(sequence[0] >= gMonoMass_x100[E] - gToleranceWide 
		&& sequence[0] <= gMonoMass_x100[E] + gToleranceWide)
	{
		nTerminalE = TRUE;
	}
	else
	{
		nTerminalE = FALSE;
	}
	
/*	Start the calculations and searches*/
	for(i = 0; i < seqLength - 1; i++)
	{
		/*Count the number of positively charged amino acids in the sequence*/
		if((sequence[i] >= gMonoMass_x100[R] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[R] + gToleranceWide) || 
			(sequence[i] >= gMonoMass_x100[H] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[H] + gToleranceWide) ||
			(sequence[i] >= gMonoMass_x100[K] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[K] + gToleranceWide))
		{
			posResidues++;
		}
		else	/*check to see if its a two aa gap that might have a positive charge*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] >= gArgPlus[j] - gToleranceWide 
					&& sequence[i] <= gArgPlus[j] + gToleranceWide) || 
					(sequence[i] >= gHisPlus[j] - gToleranceWide 
					&& sequence[i] <= gHisPlus[j] + gToleranceWide) ||
					(sequence[i] >= gLysPlus[j] - gToleranceWide 
					&& sequence[i] <= gLysPlus[j] + gToleranceWide))
				{
					posResidues++;
					break;
				}
			}
		}
		
		/*count the number of oxidized Met's, or Phe's (which could be oxidized Met)*/
		if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0  && gToleranceNarrow < 5)	/*mass accuracy sufficient
																	to determine oxMet*/
		{
			if(sequence[i] >= gMonoMass_x100[9] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[9] + gToleranceNarrow)
			{
				oxMetCount++;	/*y ions gained a oxMet*/
			}
			if(oxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}
		else	/*mass accuracy not sufficient to differentiate oxMet from Phe*/
		{
			if(sequence[i] >= gMonoMass_x100[F] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[F] + gToleranceNarrow)
			{
				oxMetCount++;	/*y ions gained a oxMet*/
			}
			if(oxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}
		
		/*Decide if this is a gap*/
		isItAGap = TRUE;

		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] <= gMonoMass_x100[j] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] - gToleranceWide)
			{
				isItAGap = FALSE;
				break;
			}
			if(sequence[i] <= gMonoMass_x100[j] + gMonoMass_x100[P] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] + gMonoMass_x100[P] - gToleranceWide)
			{
				if(proGapCount < proGapLimit)
				{
					isItAGap = FALSE;	/*if the gap could contain Pro, then don't call it a gap*/
					proGapCount++;
					break;
				}
			}
		}
		
		/*	Don't count b1 ions */
		if(!isItAGap && i == 0)
		{
			bIonTemplate += sequence[i];	/*need to add the mass to have correct b ion series later*/
			continue;
		}
		
		/*Calc b related ions assuming a single charge*/
		bIonTemplate += sequence[i];
		bIonMin17Template = bIonTemplate - ammonia;
		bIonMin18Template = bIonTemplate - water;
		bIonMin64Template = bIonTemplate - lossOf64;
		aIonTemplate = bIonTemplate - carbonMonoxide;
		for(j = 1; j <= maxCharge; j++)	/*check different charge states*/
		{
			bIon 		= (bIonTemplate + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			bIonMin17 	= (bIonMin17Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			bIonMin18 	= (bIonMin18Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			bMin64		= (bIonMin64Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			aIon		= (aIonTemplate + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			bIonTest = FALSE;
			aIonTest = FALSE;
			bMin18Test = FALSE;
			bMin17Test = FALSE;
			bMin64Test = FALSE;
			/*apply constraints to charge and mass*/
			if(bIon * j > (j-1) * 400 * gMultiplier && posResidues >= j)
			{
				/*Don't mess with the score unless b ion is less than precursor, or its LCQ data*/
				if((bIonTemplate < precursor && gParam.fragmentPattern) || gParam.fragmentPattern == 'L')
				{
					for(k = 0; k < ionCount; k++)
					{
						if(mass[k] > bIon + gToleranceWide)
						{
							break;	/*don't waste any more time looking*/
						}
						if(mass[k] > bIon - gToleranceWide)
						{
							bIonTest = TRUE;
							massDiff = fabs(mass[k] - bIon);
							bErrProb = CalcIonFound(0, massDiff);
						}
					}
					
					if(bIonTest || nTerminalQ || nTerminalE || argPresent)	/*there's a b ion, so look for the b-17 and b-18*/
					{
						for(k = 0; k < ionCount; k++)
						{
							if(bIonTest || nTerminalQ || argPresent)
							{
								if(mass[k] > bIonMin17 - gToleranceWide &&
									mass[k] < bIonMin17 + gToleranceWide)
								{
									bMin17Test = TRUE;
									massDiff = fabs(mass[k] - bIonMin17);
									b17ErrProb = CalcIonFound(0, massDiff);
								}
							}
							
							if(bIonTest || nTerminalE || argPresent)
							{
								if(mass[k] > bIonMin18 - gToleranceWide &&
									mass[k] < bIonMin18 + gToleranceWide)
								{
									bMin18Test = TRUE;
									massDiff = fabs(mass[k] - bIonMin18);
									b18ErrProb = CalcIonFound(0, massDiff);
								}
							}
							
							if(bIonTest)
							{
								if(mass[k] > aIon - gToleranceWide &&
									mass[k] < aIon + gToleranceWide)
								{
									aIonTest = TRUE;
									massDiff = fabs(mass[k] - aIon);
									aErrProb = CalcIonFound(0, massDiff);
								}
							}
							if(bIonTest && oxMetCount > 0)
							{
								if(mass[k] > bMin64 - gToleranceWide &&
									mass[k] < bMin64 + gToleranceWide)
								{
									bMin64Test = TRUE;
									massDiff = fabs(mass[k] - bMin64);
									b64ErrProb = CalcIonFound(0, massDiff);
								}
							}
						}
					}
					
					/*Calculate the probability scores*/
					/*but first find the approximate index value for the b ion (to get correct randomProb)*/
					for(k = 0; k < ionCount; k++)
					{
						if(mass[k] > bIon + gToleranceWide)
						{
							bIonIndex = k - 1;
							break;
						}
					}
					if(bIonIndex >= ionCount)	/*check the ends of the array*/
					{
						bIonIndex = ionCount - 1;
					}
					if(bIonIndex < 0)
					{
						bIonIndex = 0;
					}
					
					if(bIonTest)	/*if the calculated b ion is present*/
					{
						if(j == 1)	/*for singly charged fragments*/
						{
							individualProb = bIonProb / randomProb[bIonIndex];
							individualProb *= bErrProb;
							if(individualProb > 1)	/*ion found means normalized prob over 1 so as not to penalize*/
							{
								probScore *= individualProb;
							}
						}
						else	/*for multiply charged fragments*/
						{
							individualProb = bIonProb * bDoublyProbMultiplier / randomProb[bIonIndex];
							individualProb *= bErrProb;
							if(individualProb > 1)
							{
								probScore *= individualProb;
							}
						}
					}
					else	/*if the calculated b ion is not present*/
					{
						if(bIon > lowMassLimit && bIon < highMassLimit)	/*dont penalize if outside these limits*/
						{
							if(!(nTerminalQ && bMin17Test))	/*if b-17 is present and N-term Q, then don't penalize*/
							{
								if(!(nTerminalE && bMin18Test))	/*if b-18 is present and N-term E, then don't 
																penalize*/
								{
									if(!(argPresent && (bMin17Test || bMin18Test)))	/*if more Arg's than charges, don't penalize*/
									{
										if(j == 1)
										{
											individualProb = (1-bIonProb) / (1 - randomProb[bIonIndex]);
											if(individualProb < 1 && individualProb > 0)	/*penalize by being 
																							between 0 and 1*/
											{
												probScore *= individualProb;
											}
										}
										else
										{
											individualProb = (1 - bIonProb * bDoublyProbMultiplier) 
																/ (1 - randomProb[bIonIndex]);
											if(individualProb < 1 && individualProb > 0)
											{
												probScore *= individualProb;
											}
										}
									}
								}
							}
						}
					}
					
					/*if any of the neutral losses from b ions are present*/
					if(bMin18Test || bMin17Test || aIonTest)
					{
						if(bMin18Test)	/*if the calculated b-18 ion is present*/
						{
							if(j == 1)	/*for singly charged fragments*/
							{
								individualProb = bMinWaterProb / randomProb[bIonIndex];
								individualProb *= b18ErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}
							else	/*for multiply charged fragments*/
							{
								individualProb = bMinWaterProb * bDoublyProbMultiplier / randomProb[bIonIndex];
								individualProb *= b18ErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}		
						}
						if(bMin17Test)	/*if the calculated b-17 ion is present*/
						{
							if(j == 1)	/*for singly charged fragments*/
							{
								individualProb = bMinAmmoniaProb / randomProb[bIonIndex];
								individualProb *= b17ErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}		
							else	/*for multiply charged fragments*/
							{
								individualProb = bMinAmmoniaProb * bDoublyProbMultiplier / randomProb[bIonIndex];
								individualProb *= b17ErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}
						}
						if(aIonTest)	/*if the calculated a ion is present*/
						{
							if(j == 1)	/*for singly charged fragments*/
							{
								individualProb = aIonProb / randomProb[bIonIndex];
								individualProb *= aErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}
							else	/*for multiply charged fragments*/
							{
								individualProb = aIonProb * bDoublyProbMultiplier / randomProb[bIonIndex];
								individualProb *= aErrProb;
								if(individualProb > 1)
								{
									probScore *= individualProb;
								}
							}
						}
					}
					else	/*missing ion needs to be penalized*/
					{
						neutralLossProb = bMinWaterProb;	/*Figure out which neutral loss is least likely*/
						if(bMinAmmoniaProb < neutralLossProb)
						{
							neutralLossProb = bMinAmmoniaProb;
						}
						if(aIonProb < neutralLossProb)
						{
							neutralLossProb = aIonProb;
						}
						if(aIon > lowMassLimit && bIonMin17 < highMassLimit)
						{
							if(j == 1)
							{
								individualProb = (1 - neutralLossProb) / (1 - randomProb[bIonIndex]);
								if(individualProb < 1 && individualProb > 0)
								{
									probScore *= individualProb;
								}
							}
							else
							{
								individualProb = (1 - neutralLossProb * bDoublyProbMultiplier) / (1 - randomProb[bIonIndex]);
								if(individualProb < 1 && individualProb > 0)
								{
									probScore *= individualProb;
								}
							}
						}
					}
					if(bMin64Test)	/*don't penalize if oxMet neutral loss is absent*/
					{
						individualProb = bMin64IonProb / randomProb[bIonIndex];
						individualProb *= b64ErrProb;
						if(individualProb > 1)
						{
							probScore *= individualProb;
						}
					}
					
					if(isItAGap && j == 1 && i > 0)	/*a gap arises from lack of a b/y ion, so penalize, 
													but only do it once for j=1; also don't penalize
													a gap at the N-terminus*/
					{
						if(bIon > lowMassLimit && bIon < highMassLimit)
						{
							individualProb = (1-bIonProb) / (1 - randomProb[bIonIndex]);
							if(individualProb < 1 && individualProb > 0)
							{	
								probScore *= individualProb;
							}
						}
					}
					
					/*Now calculate the max prob score*/
					if(bIon > lowMassLimit && bIon < highMassLimit)	/*dont add to score if outside these limits*/
					{
						/*add b ion contribution*/
						if(j == 1)
						{
							individualProb = bIonProb / randomProb[bIonIndex];
							if(individualProb > 1)	
							{
								gProbScoreMax *= individualProb;
							}
						}
						else
						{
							individualProb = bIonProb * bDoublyProbMultiplier / randomProb[bIonIndex];
							if(individualProb > 1)
							{
								gProbScoreMax *= individualProb;
							}
						}
						
						/*add neutral loss contribution (only add one neutral loss per residue*/
						neutralLossProb = bMinWaterProb;	/*Figure out which neutral loss is least likely*/
						if(bMinAmmoniaProb < neutralLossProb)
						{
							neutralLossProb = bMinAmmoniaProb;
						}
						if(aIonProb < neutralLossProb)
						{
							neutralLossProb = aIonProb;
						}
						if(j == 1)	/*for singly charged fragments*/
						{
							individualProb = neutralLossProb / randomProb[bIonIndex];
							if(individualProb > 1)
							{
								gProbScoreMax *= individualProb;
							}
						}
						else	/*for multiply charged fragments*/
						{
							individualProb = neutralLossProb * bDoublyProbMultiplier / randomProb[bIonIndex];
							if(individualProb > 1)
							{
								gProbScoreMax *= individualProb;
							}
						}
						/*now add extra prob score for ions that should be present within a gap*/
						if(isItAGap && j == 1)	/*a gap arises from lack of a b/y ion, so penalize, 
													but only do it once for j=1*/
						{
							/*decide if its a legitimate 2 aa gap or something bigger*/
							TwoAAGap = FALSE;
							for(k = 0; k < gGapListIndex; k++)
							{
								if(sequence[i] <= gGapList[k] + gToleranceWide &&
									sequence[i] >= gGapList[k] - gToleranceWide)
								{
									TwoAAGap = TRUE;	/*found it as a 2aa gap*/
									break;
								}
							}
							if(TwoAAGap)
							{
								addBIons = 1;
							}
							else
							{
								addBIons = ((REAL_4)sequence[i] / (gMultiplier * AV_RESIDUE_MASS)) + 0.5;
								addBIons = addBIons - 1;	/*Ex: if three residues, then add two y ions*/
								if(addBIons < 1)
								{
									addBIons = 1;
								}
							}
							
							if(i == 0)
							{
								addBIons = addBIons - 1;	/*accounts for complete absence of b1 ions*/
							}
							
							if(bIon > lowMassLimit && bIon < highMassLimit)
							{
								individualProb = bIonProb / randomProb[bIonIndex];	
								if(individualProb > 0)
								{
									for(k = 0; k < addBIons; k++)
									{
										gProbScoreMax *= individualProb;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return(probScore);	
}

/**************************************FindYIons*********************************************
*
*	Find the y ions for the sequence and change the ionFound to 1.  Return a value that corresponds
*	to the number of consecutive y ions.
*/

REAL_8 FindYIons(INT_4 *mass, INT_4 ionCount, REAL_8 probScore,
					REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength, char argPresent)
{		
	INT_4	i, j, k, yIonIndex, posResidues, proGapLimit, proGapCount, oxMetCount, maxCharge;
	INT_4 	addYIons;
	REAL_4	water, ammonia, yIonTemplate, yIonMin17Template, yIonMin18Template;
	REAL_4	yIon, yIonMin17, yIonMin18, individualProb, neutralLossProb;
	REAL_4	massDiff, yErrProb, y17ErrProb, y18ErrProb, precursor, lowMassLimit, highMassLimit;
	REAL_4	yMin64, lossOf64, yIonMin64Template, y64ErrProb;
	BOOLEAN	yIonTest, yMin17Test, yMin18Test, isItAGap, yMin64Test, TwoAAGap;
	
	
/*	Initialize*/
	water = gElementMass_x100[OXYGEN] + gElementMass_x100[HYDROGEN] * 2;
	ammonia = gElementMass_x100[NITROGEN] + gElementMass_x100[HYDROGEN] * 3;
	lossOf64 = gElementMass_x100[CARBON] + gElementMass_x100[SULFUR] + 
				gElementMass_x100[OXYGEN] + 4 * gElementMass_x100[HYDROGEN];
	yIonTemplate = gParam.modifiedCTerm + 2 * gElementMass_x100[HYDROGEN];
	posResidues = 1;
	oxMetCount = 0;
	if(gParam.peptideMW < 1000 * gMultiplier)
	{
		proGapLimit = 1;
	}
	else
	{
		proGapLimit = 2;
	}
	proGapCount = 0;
	if(gParam.chargeState == 1)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = 2;
	}
	
	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;	
	if(gParam.fragmentPattern == 'L')	/*ion masses outside of this range are not penalized if not found*/
	{
		lowMassLimit = precursor * 0.333;	/*so-called 1/3 rule*/
		highMassLimit = 2000 * gMultiplier;	/*mass limit for Deca*/
	}
	else
	{
		lowMassLimit = 146 * gMultiplier;	/*y1 for Lys*/
		highMassLimit = 2 * precursor;	/*often the very high mass ions are missing*/
	}
	
	for(i = seqLength - 1; i > 0; i--)
	{
		/*Determine if its a positive residue*/
		if((sequence[i] >= gMonoMass_x100[R] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[R] + gToleranceWide) || 
			(sequence[i] >= gMonoMass_x100[H] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[H] + gToleranceWide) ||
			(sequence[i] >= gMonoMass_x100[K] - gToleranceWide 
			&& sequence[i] <= gMonoMass_x100[K] + gToleranceWide))
		{
			posResidues++;
		}
		else	/*check to see if its a two aa gap that might have a positive charge*/
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if((sequence[i] >= gArgPlus[j] - gToleranceWide 
					&& sequence[i] <= gArgPlus[j] + gToleranceWide) || 
					(sequence[i] >= gHisPlus[j] - gToleranceWide 
					&& sequence[i] <= gHisPlus[j] + gToleranceWide) ||
					(sequence[i] >= gLysPlus[j] - gToleranceWide 
					&& sequence[i] <= gLysPlus[j] + gToleranceWide))
				{
					posResidues++;
					break;
				}
			}
		}
		
		/*Determine if this is a gap*/
		isItAGap = TRUE;

		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(sequence[i] < gMonoMass_x100[j] + gToleranceWide &&
				sequence[i] > gMonoMass_x100[j] - gToleranceWide)
			{
				isItAGap = FALSE;
				break;
			}
			if(sequence[i] <= gMonoMass_x100[j] + gMonoMass_x100[P] + gToleranceWide &&
				sequence[i] >= gMonoMass_x100[j] + gMonoMass_x100[P] - gToleranceWide)
			{
				if(proGapCount < proGapLimit)
				{
					isItAGap = FALSE;	/*if the gap could contain Pro, then don't call it a gap*/
					proGapCount++;
					break;
				}
			}
		}
		
		/*count the number of oxidized Met's, or Phe's (which could be oxidized Met)*/
		if(gParam.fragmentPattern == 'Q' && gParam.qtofErr != 0  && gToleranceNarrow < 5)	/*mass accuracy sufficient
																	to determine oxMet*/
		{
			if(sequence[i] >= gMonoMass_x100[9] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[9] + gToleranceNarrow)
			{
				oxMetCount++;	/*y ions gained a oxMet*/
			}
			if(oxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}
		else	/*mass accuracy not sufficient to differentiate oxMet from Phe*/
		{
			if(sequence[i] >= gMonoMass_x100[F] - gToleranceNarrow
				&& sequence[i] <= gMonoMass_x100[F] + gToleranceNarrow)
			{
				oxMetCount++;	/*y ions gained a oxMet*/
			}
			if(oxMetCount < 0)
			{
				printf("LutefiskScore:FindABYIons The number of oxidized Mets went negative.");
				exit(1);
			}
		}
		
		yIonTemplate += sequence[i];
		yIonMin17Template = yIonTemplate - ammonia;
		yIonMin18Template = yIonTemplate - water;
		yIonMin64Template = yIonTemplate - lossOf64;
		for(j = 1; j <= maxCharge; j++)	/*check different charge states*/
		{
			
			yIon 		= (yIonTemplate + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			yIonMin17 	= (yIonMin17Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			yIonMin18 	= (yIonMin18Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			yMin64		= (yIonMin64Template + (j-1)*gElementMass_x100[HYDROGEN]) / j;
			yIonTest = FALSE;
			yMin17Test = FALSE;
			yMin18Test = FALSE;
			yMin64Test = FALSE;
			
			/*apply constraints to charge and mass*/
			if(yIon * j > (j-1) * 400 * gMultiplier && posResidues >= j)
			{
				for(k = 0; k < ionCount; k++)
				{
					if(mass[k] > yIon + gToleranceWide)
					{
						break;	/*don't waste any more time looking*/
					}
					if(mass[k] > yIon - gToleranceWide)
					{
						yIonTest = TRUE;
						massDiff = fabs(mass[k] - yIon);
						yErrProb = CalcIonFound(0, massDiff);
					}
				}
				if(yIonTest || argPresent)	/*there's a y ion, so look for the y-17 and y-18, or for non-mobile
											proton fragmentation*/
				{
					for(k = 0; k < ionCount; k++)
					{
						if(mass[k] > yIonMin17 - gToleranceWide &&
							mass[k] < yIonMin17 + gToleranceWide)
						{	/*y-17 intensity should be less than the y ion*/
							yMin17Test = TRUE;
							massDiff = fabs(mass[k] - yIonMin17);
							y17ErrProb = CalcIonFound(0, massDiff);
						}
						if(mass[k] > yIonMin18 - gToleranceWide &&
							mass[k] < yIonMin18 + gToleranceWide)
						{
							yMin18Test = TRUE;
							massDiff = fabs(mass[k] - yIonMin18);
							y18ErrProb = CalcIonFound(0, massDiff);
						}
						if(oxMetCount > 0)
						{
							if(mass[k] > yMin64 - gToleranceWide &&
								mass[k] < yMin64 + gToleranceWide)
							{
								yMin64Test = TRUE;
								massDiff = fabs(mass[k] - yMin64);
								y64ErrProb = CalcIonFound(0, massDiff);
							}
						}
					}
				}
				
				
				/*Calculate the probability scores*/
				/*but first find the approximate index value for the y ion (to get the correct randomProb)*/
				for(k = 0; k < ionCount;k++)
				{
					if(mass[k] < yIon + gToleranceWide)
					{
						yIonIndex = k - 1;
						break;
					}
				}
				if(yIonIndex >= ionCount)
				{
					yIonIndex = ionCount - 1;
				}
				if(yIonIndex < 0)
				{
					yIonIndex = 0;
				}
				
				if(yIonTest)	/*if the calculated y ion is present*/
				{
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = yIonProb / randomProb[yIonIndex];
						individualProb *= yErrProb;
						if(individualProb > 1)
						{
							probScore *= individualProb;
						}
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = yIonProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						individualProb *= yErrProb;
						if(individualProb > 1)
						{
							probScore *= individualProb;
						}
					}
				}
				else	/*if the calculated y ion is not present*/
				{
					if(yIon > lowMassLimit && yIon < highMassLimit && i > 1)
					{
						if(!(argPresent && (yMin17Test || yMin18Test)))
						{
							if(j == 1)
							{
								individualProb = (1-yIonProb) / (1 - randomProb[yIonIndex]);
								if(individualProb < 1 && individualProb > 0)
								{
									probScore *= individualProb;
								}
							}
							else
							{
								individualProb = (1 - yIonProb * yDoublyProbMultiplier) / (1 - randomProb[yIonIndex]);
								if(individualProb < 1 && individualProb > 0)
								{
									probScore *= individualProb;
								}
							}
						}
					}
				}

				if(yMin18Test || yMin17Test)	/*if the calculated y-18 ion is present*/
				{
					if(yMin18Test)
					{
						if(j == 1)	/*for singly charged fragments*/
						{
							individualProb = yMinWaterProb / randomProb[yIonIndex];
							individualProb *= y18ErrProb;
							if(individualProb > 1)
							{
								probScore *= individualProb;
							}
						}
						else	/*for multiply charged fragments*/
						{
							individualProb = yMinWaterProb * yDoublyProbMultiplier / randomProb[yIonIndex];
							individualProb *= y18ErrProb;
							if(individualProb > 1)
							{
								probScore *= individualProb;
							}
						}
					}
					if(yMin17Test)	/*if the calculated y-18 ion is present*/
					{
						if(j == 1)	/*for singly charged fragments*/
						{
							individualProb = yMinAmmoniaProb / randomProb[yIonIndex];
							individualProb *= y17ErrProb;
							if(individualProb > 1)
							{
								probScore *= individualProb;
							}
						}
						else	/*for multiply charged fragments*/
						{
							individualProb = yMinAmmoniaProb * yDoublyProbMultiplier / randomProb[yIonIndex];
							individualProb *= y17ErrProb;
							if(individualProb > 1)
							{
								probScore *= individualProb;
							}
						}
					}
				}
				else	/*missing ions need to be penalized*/
				{
					if(yMinWaterProb < yMinAmmoniaProb)
					{
						neutralLossProb = yMinWaterProb;
					}
					else
					{
						neutralLossProb = yMinAmmoniaProb;
					}
					if(yIonMin18 > lowMassLimit && yIonMin17 < highMassLimit && i > 1)
					{
						if(j == 1)
						{
							individualProb = (1 - neutralLossProb) / (1 - randomProb[yIonIndex]);
							if(individualProb < 1 && individualProb > 0)
							{
								probScore *= individualProb;
							}
						}
						else
						{
							individualProb = (1 - neutralLossProb * yDoublyProbMultiplier) / (1 - randomProb[yIonIndex]);
							if(individualProb < 1 && individualProb > 0)
							{
								probScore *= individualProb;
							}
						}
					}
				}
				if(yMin64Test)	/*don't penalize if oxMet neutral loss is absent*/
				{
					individualProb = yMin64IonProb / randomProb[yIonIndex];
					individualProb *= y64ErrProb;
					if(individualProb > 1)
					{
						probScore *= individualProb;
					}
				}
				
				if(isItAGap && j == 1 && i > 0)	/*a gap arises from lack of a y ion, so penalize, 
															but only do it once for j=1; also, don't penalize
															for a gap at the N-terminus*/
				{
					if(yIon > lowMassLimit && yIon < highMassLimit)
					{
						individualProb = (1-yIonProb) / (1 - randomProb[yIonIndex]);	
						if(individualProb < 1 && individualProb > 0)
						{
							probScore *= individualProb;
						}
					}
				}
				/*calculate y ion contribution to max possible prob score for this sequence*/
				if(yIon > lowMassLimit && yIon < highMassLimit)
				{
					/*Add y ion contribution*/
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = yIonProb / randomProb[yIonIndex];
						if(individualProb > 1)
						{
							gProbScoreMax *= individualProb;
						}
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = yIonProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						if(individualProb > 1)
						{
							gProbScoreMax *= individualProb;
						}
					}
					/*add one neutral loss contribution*/
					if(yMinWaterProb < yMinAmmoniaProb)
					{
						neutralLossProb = yMinWaterProb;
					}
					else
					{
						neutralLossProb = yMinAmmoniaProb;
					}
					if(j == 1)	/*for singly charged fragments*/
					{
						individualProb = neutralLossProb / randomProb[yIonIndex];
						if(individualProb > 1)
						{
							gProbScoreMax *= individualProb;
						}
					}
					else	/*for multiply charged fragments*/
					{
						individualProb = neutralLossProb * yDoublyProbMultiplier / randomProb[yIonIndex];
						if(individualProb > 1)
						{
							gProbScoreMax *= individualProb;
						}
					}
					if(isItAGap && j == 1)	/*a gap arises from lack of a y ion, so penalize, 
															but only do it once for j=1*/
					{
						/*decide if its a legitimate 2 aa gap or something bigger*/
						TwoAAGap = FALSE;
						for(k = 0; k < gGapListIndex; k++)
						{
							if(sequence[i] <= gGapList[k] + gToleranceWide &&
								sequence[i] >= gGapList[k] - gToleranceWide)
							{
								TwoAAGap = TRUE;	/*found it as a 2aa gap*/
								break;
							}
						}
						if(TwoAAGap)
						{
							addYIons = 1;
						}
						else
						{
							addYIons = ((REAL_4)sequence[i] / (gMultiplier * AV_RESIDUE_MASS)) + 0.5;
							addYIons = addYIons - 1;	/*Ex: if three residues, then add two y ions*/
							if(addYIons < 1)
							{
								addYIons = 1;
							}
						}
						if(yIon > lowMassLimit && yIon < highMassLimit)
						{
							individualProb = yIonProb / randomProb[yIonIndex];	
							if(individualProb > 0)
							{
								for(k = 0; k < addYIons; k++)
								{
									gProbScoreMax *= individualProb;
								}
							}
						}
					}
				}
			}
		}
	}

	/*Check to see if the N-terminal residue is a gap greater than 2aa's*/
	TwoAAGap = FALSE;
	for(k = 0; k < gGapListIndex; k++)
	{
		if(sequence[0] <= gGapList[k] + gToleranceWide &&
			sequence[0] >= gGapList[k] - gToleranceWide)
		{
			TwoAAGap = TRUE;	/*found it as a 2aa gap or 1aa residue*/
			break;
		}
	}
	if(!TwoAAGap)
	{

		addYIons = ((REAL_4)sequence[0] / (gMultiplier * AV_RESIDUE_MASS)) + 0.5;
		addYIons = addYIons - 1;	/*Ex: if three residues, then add two y ions*/
		if(addYIons < 1)
		{
			addYIons = 1;
		}
		if(addYIons > 0)
		{
			individualProb = yIonProb / randomProb[yIonIndex];	
			if(individualProb > 0)
			{
				for(k = 0; k < addYIons; k++)
				{
					gProbScoreMax *= individualProb;
				}
			}
		}
	}

	return(probScore);
}




/****************************PrintScoreDetailsToXLFile***************************************************
*	This function prints header information to the output file.
*/

void PrintScoreDetailsToXLFile(struct SequenceScore *firstScorePtr, REAL_4 perfectProbScore)
{
	FILE *fp;
	INT_4 i, j, seqNum;
	REAL_4 xcorrNormalizer;
	struct SequenceScore *maxPtr, *currPtr;
   	const 	time_t		theTime = (const time_t)time(NULL);
	char  outputFile[256], fileName[256];
    INT_4 length;
    INT_4 fileCount;
    
	/*Make up a name for the file*/
 	if (strlen(gParam.cidFilename) != 0)
    {
       

        /* Start from the CID filename */
        strcpy (outputFile, gParam.cidFilename);

        length = strlen(outputFile);

        strcat(outputFile, ".xl");

    }
    else
    {
    	printf("not printing details");
    	return;
    }
    
 /* Make sure that the file doesn't already exist. If it does, append a number. */
 
 
	strcpy(fileName, outputFile);
	fileCount = 1;

	while (1)
	{
		FILE *fp = fopen(fileName, "r");

		if (NULL == fp) break;

		fclose(fp);

		strcpy(fileName, outputFile);
		sprintf(fileName + strlen(fileName), "%d\0", fileCount++);

		if (fileCount > 20)
		{
			printf("Too many old output files! Please clean up a bit first! Quitting.");
			exit(1);
		}
	}     
        
        
    /* Open a new file.*/
	fp = fopen(fileName, "w");
	if(fp == NULL)	/*fopen returns NULL if there's a problem.*/
	{
		printf("Cannot open %s to write the output.\n", gParam.outputFile);
		exit(1);
	}

	fprintf(fp, "Run Date: %20s", ctime(&theTime));

        /* Print header information from gParam to the console and the file.*/
	fprintf(fp, " Filename: ");	/*Print the CID data file name.*/
	i = 0;
	while(gParam.cidFilename[i] != 0)
	{
		fputc(gParam.cidFilename[i], fp);
		i++;
	}
	fprintf(fp, "\n Molecular Weight: %7.2f", gParam.peptideMW);
	fprintf(fp, "  Molecular Weight Tolerance: %5.2f", gParam.peptideErr);
	fprintf(fp, "  Fragment Ion Tolerance: %5.2f", gParam.fragmentErr);
	fprintf(fp, "\n Ion Offset: %5.2f", gParam.ionOffset);
	fprintf(fp, "  Charge State: %2ld", gParam.chargeState);
	if(gParam.centroidOrProfile == 'P')
	{
		fprintf(fp, "   Profile Data \n");
	}
	else
	{
		fprintf(fp, "   Centroided or Pre-processed Data \n");
	}
	if(gParam.proteolysis == 'T')
	{
		fprintf(fp, " Tryptic Digest");
	}
	if(gParam.proteolysis == 'K')
	{
		fprintf(fp, " Lys-C Digest");
	}
	if(gParam.proteolysis == 'E')
	{
		fprintf(fp, " Glu-C Digest");
	}
	if(gParam.proteolysis == 'N')
	{
		fprintf(fp, " ??? Digest");
	}
	if(gParam.fragmentPattern == 'G')
	{
		fprintf(fp, "       Unknown Fragmentation Pattern \n");
	}
	if(gParam.fragmentPattern == 'T')
	{
		fprintf(fp, "       Tryptic Triple Quadrupole Fragmentation Pattern \n");
	}
	if(gParam.fragmentPattern == 'L')
	{
		fprintf(fp, "       Tryptic Ion Trap Fragmentation Pattern \n");
	}
	if(gParam.fragmentPattern == 'Q')
	{
		fprintf(fp, "       Tryptic QTOF Fragmentation Pattern \n");
	}
	fprintf(fp, " Cysteine residue mass: %7.2f", gParam.cysMW);
	fprintf(fp, "  Switch from monoisotopic to average mass at %d \n", gParam.monoToAv);
	fprintf(fp, " Ions per window: %.1f", gParam.ionsPerWindow);
	fprintf(fp, "  Extension Threshold: %4.2f", gParam.extThresh);
	fprintf(fp, "  Extension Number: %2ld", gParam.maxExtNum);
	fprintf(fp, "\n Gaps: %2ld", gParam.maxGapNum);
	fprintf(fp, " Peak Width: %4.1f", ((gParam.peakWidth) * 2));
	fprintf(fp, " Data Threshold: %5.2f (%ld)", gParam.ionThreshold, gParam.intThreshold);
	fprintf(fp, " Ions per residue: %.1f", gParam.ionsPerResidue);
	fprintf(fp, "\n Amino acids known to be present: ");
	i = 0;
	while(gParam.aaPresent[i] != 0)
	{
		fputc(gParam.aaPresent[i], fp);
		i++;
	}
	fprintf(fp, "\n Amino acids known to be absent: ");
	i = 0;
	while(gParam.aaAbsent[i] != 0)
	{
		fputc(gParam.aaAbsent[i], fp);
		i++;
	}
	fprintf(fp, "\n");
	fprintf(fp, "\n C-terminal mass: %7.4f", gParam.modifiedCTerm);
	fprintf(fp, "\n N-terminal mass: %7.4f", gParam.modifiedNTerm);
	fprintf(fp, "\n N-terminal Tag Mass: %7.2f", gParam.tagNMass);
	fprintf(fp, " C-terminal Tag Mass: %7.2f", gParam.tagCMass);
	fprintf(fp, " Sequence Tag: ");
	i = 0;
	while(gParam.tagSequence[i] != 0)
	{
		fputc(gParam.tagSequence[i], fp);
		i++;
	}
	if(gParam.edmanPresent)
	{
		fprintf(fp, "\n Edman data is available.    ");
	}
	else
	{
		fprintf(fp, "\n Edman data is not available.    ");
	}
	if(gParam.autoTag == TRUE)
	{
		fprintf(fp, "AutoTag ON");
	}
	else
	{
		fprintf(fp, "AutoTag OFF");
	}
	if(gParam.CIDfileType == 'T')
	{
		fprintf(fp, "   CID data file is tab-delineated");
	}
	if(gParam.CIDfileType == 'F')
	{
		fprintf(fp, "   CID data file is Finnigan ASCII file");
	}
	 
	/*Print the details to the xl file*/
	fprintf(fp, "\n Perfect probability score: %6.2f", perfectProbScore);
	
/*	Find the xcorr normalizer.*/
	xcorrNormalizer = 0;
	currPtr = firstScorePtr;
	while(currPtr != NULL)
	{
		if(currPtr->crossDressingScore > xcorrNormalizer)
		{
			xcorrNormalizer = currPtr->crossDressingScore;
		}
		currPtr = currPtr->next;
	}
	xcorrNormalizer = 1;
/*	Set up the screen to print some of the output.*/
	fprintf(fp, "\n Rank  X-corr  IntScr  IntOnlyScr Quality  ProbScr StDevErr CS   CalFact  Sequence\n");
	
/*	Count the sequences.*/
	seqNum = 0;
	maxPtr = firstScorePtr;
	while(maxPtr != NULL)
	{
		seqNum++;
		maxPtr = maxPtr->next;
	}
		
	for(i = 1; i <= 500 && i <= seqNum; i++)	/*List the top 500 sequences.*/
	{
		maxPtr = firstScorePtr;
		while(maxPtr != NULL)
		{	
			if(maxPtr->rank == i)
			{
				/*Change peptide[j] to single letter code.*/
				char *peptideString;
				INT_4 peptide[MAX_PEPTIDE_LENGTH];
				INT_4 peptideLength = 0;
				
				j = 0;
				while(maxPtr->peptide[j] != 0)
				{
					peptide[j] = maxPtr->peptideSequence[j];
					peptideLength++;
					j++;
				}
				peptideString = PeptideString(peptide, peptideLength);
				if(maxPtr->databaseSeq)
				{
					strcat(peptideString, " ");	/*used to denote this was database sequence*/
				}
				if(peptideString) 
				{
					fprintf(fp, " %3ld   %5.3f   %5.3f   %5.3f      %5.3f    %5.3f  %6.4f  %2ld   %8.6f %s\n", i,
						 maxPtr->crossDressingScore / xcorrNormalizer, maxPtr->intensityScore,
						 maxPtr->intensityOnlyScore, maxPtr->quality, maxPtr->probScore, 
						 maxPtr->stDevErr, maxPtr->cleavageSites, 
						 maxPtr->calFactor, peptideString);
					free(peptideString);
				}
			
				break;
			}
			
			maxPtr = maxPtr->next;
		}
	}


	fclose(fp);
	return;
}


/***********************************CalcPerfectProbScore******************************************
*
*	Determine a theoretical perfect probScore, given the list of ions present.
*/
REAL_4	CalcPerfectProbScore(INT_4 fragNum, INT_4 *fragMOverZ)
{
	REAL_8 score = 0.95; 
	REAL_4 individualProb, probScore;
	REAL_4 massRange = fragMOverZ[fragNum-1] - fragMOverZ[0];
	REAL_4 randomMatchProb = (REAL_4)fragNum / (REAL_4)massRange * gMultiplier;
	REAL_4 residueMass = AV_RESIDUE_MASS * gMultiplier;
	INT_4 seqLength = (gParam.peptideMW * gMultiplier / residueMass) + 0.5;
	INT_4 halfSeqLength = seqLength / 2;
	INT_4 i;
	
	for(i = seqLength - 1; i > 1; i--)
	{
		individualProb = yIonProb / randomMatchProb;
		score *= individualProb;
	}
	for(i = 2; i < seqLength - 1; i++)
	{
		if(i < halfSeqLength || gParam.fragmentPattern =='L')
		{
			individualProb = bIonProb / randomMatchProb;
			score *= individualProb;
		}
	}
	
	if(score > 1)
	{
		probScore = log10(score);
	}
	else	/*keep things positive by only logging things over a value of 1*/
	{
		probScore = 0;
	}
	
	/*normalize by length*/
	probScore = probScore / (2 * seqLength);
	
	return(probScore);
}
