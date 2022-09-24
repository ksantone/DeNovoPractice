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



/* ANSI headers */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* Haggis headers */
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"
#include "ListRoutines.h"

#if(defined(__MWERKS__) && __dest_os == __mac_os)
#include "getopt.h"

#include <Types.h>
#include <StandardFile.h>		
#include <sioux.h>
#include <console.h>			
	StandardFileReply freply;	
	Point wpos;			
	INT_4 tval;			
	char prompt[256];		
#endif

#if(defined(__MWERKS__) && __dest_os == __win32_os)
#include "getopt.h"
#endif

/*Definitions for this file*/
//#define MIN_NUM_IONS 			5		/*Minimum number of ions after processing in GetCID*/
//#define MAX_ION_MASS 			3000	/*Ions greater than this are deemed too high to not be a mistake*/
//#define MIN_HIGHMASS_INT_RATIO 	0.1		/*Ratio of high mass intensity over total intensity*/
//#define HIGH_MASS_RATIO 		0.9		/*Ions are counted until this % of high mass ion intensity is reached*/
//#define LCQ_INT_SUM_CUTOFF 		500		/*Cutoff for good intensity total for LCQ data*/
//#define QTOF_INT_SUM_CUTOFF		140  	/*Cutoff for good intensity total for Qtof data*/
//#define MAX_HIGH_MASS 			100		/*Max number of ions greater than precursor*/
//#define MAX_ION_NUM				200		/*Max number of ions*/
#define MAX_SEQUENCES	10000	/*Max number of sequences to store*/
//#define MAX_MASS 				2500	/*Peptides above this mass are tossed out.*/
//#define MIN_MASS 				800		/*Peptides below this mass are tossed out.*/
//#define LOW_MASS_ION_NUM 		19		/*Number of peptide-related low mass ions*/

/*Global variables for this file*/
INT_4 gForwardNodeConnect[MAX_ION_NUM][AMINO_ACID_NUMBER];
INT_4 gBackwardNodeConnect[MAX_ION_NUM][AMINO_ACID_NUMBER];
INT_4 gForwardNum[MAX_ION_NUM], gBackwardNum[MAX_ION_NUM];
INT_4 gIonCount;
INT_4 gEdgeNum;
INT_4 gSequenceNodes[MAX_SEQUENCES][MAX_ION_NUM];
INT_4 gSeqCount;
INT_4 gSequenceNum;
INT_4 gPepLength[MAX_SEQUENCES*2];
INT_4 gPepMassSeq[MAX_SEQUENCES*2][MAX_PEPTIDE_LENGTH];
INT_4 gMatchSeries[MAX_SEQUENCES*2];
INT_4 gAAArray[AMINO_ACID_NUMBER];
INT_4 gAAMonoArray[AMINO_ACID_NUMBER];
INT_4 gAANum, gMassRange, gCTermKIndex, gCTermRIndex, gLutefiskSequenceCount;
BOOLEAN gNotTooManySequences = TRUE;

/*****************************Haggis*************************************************
*
*	First, Haggis decides how many fragment ion charge states to consider (+1 or +2,
*		rejecting precursors over +3).
*	Second, For each charge state Haggis converts the linked list firstMassPtr to a ion mass array
*		of singly-charged fragments.
*	Third, for each charge state it finds all singly-charged ions that can be connected to each other
*		via single amino acid residue mass jumps.
*	Fourth, it converts the sequences of nodes to sequences of residue masses, assuming that each
*		sequence of nodes could be either b or y ions.
*	Fifth, it adds the sequences to the linked list of sequences already produced by subsequencing.
*
*/

struct Sequence *Haggis(struct Sequence *firstSequencePtr , struct MSData *firstMassPtr)
{
	INT_4	*mass;
	INT_4	j, i, maxCharge;
	INT_4 peptide[MAX_PEPTIDE_LENGTH];
	INT_4 peptideLength;
	INT_4 score;
	INT_4 nodeValue;
	INT_2 nodeCorrection;
	INT_4 gapNum, lutefiskSequenceCount;
	struct Sequence *currPtr;
	
	/*Count and report the number of Lutefisk-derived sequences*/
	lutefiskSequenceCount = 0;
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		lutefiskSequenceCount++;
		currPtr = currPtr->next;
	}
	printf("Lutefisk sequences: %ld \n", lutefiskSequenceCount);
	gLutefiskSequenceCount = lutefiskSequenceCount;	/*need to be global for StoreSeq*/
	
	/*Don't bother working on precursor charge states more than 3*/
	if(gParam.chargeState > 3)
	{
		return(firstSequencePtr);
	}
	
	/*Determine maximum charge state of fragment ions.  Precursors of +3 have max charge of 2, +1 and +2 
	have a max charge of only 1*/
	if(gParam.chargeState == 3)
	{
		maxCharge = 2;
	}
	else
	{
		maxCharge = 1;
	}
	
	/*	Make some space*/
	mass = malloc(MAX_ION_NUM * sizeof(REAL_4));
	if(mass == NULL)
	{
		printf("Haggis:  Out of memory");
		exit(1);
	}
	
	/*Initialize variables*/
	gSequenceNum = 0;
	for(i = 0; i < MAX_SEQUENCES * 2; i++)
	{
		gPepLength[i] = 0;
		for(j = 0; j < MAX_PEPTIDE_LENGTH; j++)
		{
			gPepMassSeq[i][j] = 0;
		}
	}
	
	
	/*Consider different charge states for fragment ions*/
	for(j = 1; j <= maxCharge; j++)
	{
		/*Load mass arrays*/
		mass = LoadMassArrays(mass, firstMassPtr, j);
	
		/*Set up the backward and forward node connections*/
		SetupBackwardAndForwardNodes(mass);
		
		/*Find sequences of nodes*/
		FindNodeSequences(mass);
				
		/*Convert sequences of nodes to sequences of residue masses assuming they are
		both b and y ions*/
		GetSequenceOfResidues(mass);
	}
	
/*	Try to connect sequences.*/

	AppendSequences();
	
/*	Try to fill in the unsequenced ends with reasonable sequences. */
	
	FleshOutSequenceEnds(firstMassPtr);
	
/*	To be consistent with the Lutefisk sequences, replace sequence regions that are unsupported by y/b ions w/ 
	bracketed masses.*/
	
	ModifyHaggisSequences(firstMassPtr);
	

	/*Find the highest score in the linked list*/
	score = 0;
	currPtr = firstSequencePtr;
	while(currPtr != NULL)
	{
		if(currPtr->score > score)
		{
			score = currPtr->score;
		}
		currPtr = currPtr->next;
	}
	if(score == 0)
	{
		score = 1;	/*if all the Lutefisk sequences have score of zero, then give Haggis sequences a non-zero score*/
	}
	
	/*Assign some values for the linked list*/
	nodeValue = gParam.peptideMW - gParam.modifiedCTerm + 0.5;
	nodeCorrection = 0;
	gapNum = 0;
	
	/*Add sequences to linked list*/
	for(i = 0; i < gSequenceNum; i++)
	{
		peptideLength = gPepLength[i];
		
		if(peptideLength < MAX_PEPTIDE_LENGTH && peptideLength > 3)	/*toss out anything too small or big*/
		{
			for(j= 0; j < peptideLength; j++)
			{
				peptide[j] = gPepMassSeq[i][j];
			}
			
			firstSequencePtr = LinkHaggisSubsequenceList(firstSequencePtr, LoadHaggisSequenceStruct(peptide, 
								peptideLength, score, nodeValue, gapNum, nodeCorrection));
		}
	}
	
	printf("Haggis sequences: %ld \n", gSequenceNum);
	
	free(mass);
	return(firstSequencePtr);
}

/********************************ModifyHaggisSequences**********************************************
*
*	Search each sequence to see if there is a b or y ion between each amino acid.  If not, combine
*	the amino acids for which no evidence is available.
*/

void	ModifyHaggisSequences(struct MSData *firstMassPtr)
{
	INT_4 i, j, k;
	BOOLEAN bIonTest, yIonTest;
	
	for(i = 0; i < gSequenceNum; i++)
	{
		for(j = 0; j < gPepLength[i] - 1; j++)
		{
			bIonTest = FindBIon(i,j, firstMassPtr);
			yIonTest = FindYIon(i,j, firstMassPtr);
			if(!bIonTest && !yIonTest)
			{
				gPepMassSeq[i][j] += gPepMassSeq[i][j + 1];
				for(k = j + 1; k < gPepLength[i]; k++)
				{
					gPepMassSeq[i][k] = gPepMassSeq[i][k + 1];
				}
				gPepLength[i] -= 1;
				j--;
			}
		}
	}
	
	return;
}

/*********************************FindBIon***********************************************************
*
*	Look for a b ion.  If found, return a TRUE value.
*/

BOOLEAN	FindBIon(INT_4 sequenceIndex,INT_4 residueIndex, struct MSData *firstMassPtr)
{
	INT_4 	i;
	INT_4	bIon, maxCharge, bIon1Charge;
	BOOLEAN bIonPresent;
	struct MSData *currPtr;
	
	
	/*Initialize*/
	bIonPresent = FALSE;
	bIon1Charge = gParam.modifiedNTerm;
	if(gParam.chargeState > 2)
	{
		maxCharge = 2;
	}
	else
	{
		maxCharge = 1;
	}
	
	/*Calculate singly-charged b ion mass*/
	for(i = 0; i <= residueIndex; i++)
	{
		bIon1Charge += gPepMassSeq[sequenceIndex][i];
	}
	
	/*For each charge, look for the b ion*/
	for(i = 1; i <= maxCharge; i++)
	{
		bIon = (bIon1Charge + (i - 1) * gElementMass_x100[HYDROGEN]) / i;
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= bIon - gParam.fragmentErr * 1.5)	/*tolerance is very wide here*/
			{
				if(currPtr->mOverZ > bIon + gParam.fragmentErr * 1.5)
				{
					break;
				}
				bIonPresent = TRUE;
			}
			currPtr = currPtr->next;
		}
	}
	

	return(bIonPresent);
}

/*********************************FindYIon***********************************************************
*
*	Look for a y ion.  If found, return a TRUE value.
*/

BOOLEAN	FindYIon(INT_4 sequenceIndex,INT_4 residueIndex, struct MSData *firstMassPtr)
{
	INT_4 	i;
	INT_4	yIon, maxCharge, yIon1Charge;
	BOOLEAN yIonPresent;
	struct MSData *currPtr;
	
	
	/*Initialize*/
	yIonPresent = FALSE;
	yIon1Charge = gParam.modifiedCTerm + 2 * gElementMass_x100[HYDROGEN];
	if(gParam.chargeState > 2)
	{
		maxCharge = 2;
	}
	else
	{
		maxCharge = 1;
	}
	
	/*Calculate singly-charged y ion mass*/
	for(i = gPepLength[sequenceIndex] - 1; i > residueIndex; i--)
	{
		yIon1Charge += gPepMassSeq[sequenceIndex][i];
	}
	
	/*For each charge, look for the y ion*/
	for(i = 1; i <= maxCharge; i++)
	{
		yIon = (yIon1Charge + (i - 1) * gElementMass_x100[HYDROGEN]) / i;
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= yIon - gParam.fragmentErr * 1.5)	/*tolerance is very wide here*/
			{
				if(currPtr->mOverZ > yIon + gParam.fragmentErr * 1.5)
				{
					break;
				}
				yIonPresent = TRUE;
			}
			currPtr = currPtr->next;
		}
	}
	

	return(yIonPresent);
}

/*********************************AppendSequences****************************************************
*
*	This function finds two sequences that could be combined into a new one.  
*/

void	AppendSequences()
{
	INT_4	i, j, k, testMass, massDiff, newSeqNum, massLimit;
	INT_4	newSeqIndex[MAX_SEQUENCES][2];
	REAL_8 	maxSequences = MAX_SEQUENCES;
	BOOLEAN	storeIt, saveAll;
	
	massLimit = gMonoMass_x100[A] * 2 - gParam.fragmentErr;	/*unsequenced mass separating the N- and C-terminal
															sequences has to be more than this value*/
	
	/*Decide if all appended sequences can be saved*/
	maxSequences = 	sqrt(maxSequences);													
	if(maxSequences > gSequenceNum)
	{
		saveAll = TRUE;	/*there are too few sequences to have to worry about generating too many new ones
						so go ahead and keep all of them*/
	}
	else
	{
		saveAll = FALSE;	/*need to save only the ones that correspond to certain masse*/
	}
	
	/*Start the search for new sequences derived by sticking two old ones together*/
	newSeqNum = 0;
	for(i = 0; i < gSequenceNum; i++)
	{
		testMass = 0;
		for(j = 0; j < gPepLength[i] - 1; j++)
		{
			testMass += gPepMassSeq[i][j];	/*Find N-terminal mass plus sequence of first one*/
		}
		for(j = i; j < gSequenceNum; j++)
		{
			massDiff = gPepMassSeq[j][0] - testMass; /*Subtract C-terminal unsequenced mass from
														N-terminal mass of first one */
			if(massDiff > massLimit)	/*At moment only requirement is that mass diff be more than 2xAla*/
			{
				if(gMatchSeries[j] != i && newSeqNum < MAX_SEQUENCES)	/*don't append a sequence that is the 
																		reverse of itself*/
				{
					storeIt = FALSE;	/*assume the worst*/
					if(saveAll)
					{
						storeIt = TRUE;	/*not enough sequences to worry about overflow*/
					}
					else
					{
			//			if(massDiff > gMonoMass_x100[W] * 2 && massDiff < 400 * gMultiplier)
			//			{
			//				storeIt = FALSE;	/*if the mass diff is between 372 and 500, save it*/
			//			}
			//			else
			//			{
							for(k = 0; k < gGapListIndex; k++)
							{
								if(massDiff <= gGapList[k] + gParam.fragmentErr &&
									massDiff >= gGapList[k] - gParam.fragmentErr)
								{
									storeIt = TRUE;	/*if its a one or two amino acid mass, save it*/
									break;
								}
							}
			//			}
					}
					if(storeIt)	/*Store the index values of the two old sequences*/
					{
						newSeqIndex[newSeqNum][0] = i;	/*N-terminal bit*/
						newSeqIndex[newSeqNum][1] = j;	/*C-terminal bit*/
						newSeqNum++;
					}
				}
			}
			
		}
	}
	
	/*Make the new sequences and add them to the global list*/
	for(i = 0; i < newSeqNum; i++)
	{
		if(gSequenceNum < MAX_SEQUENCES * 2)
		{
			testMass = 0;
			for(j = 0; j < gPepLength[newSeqIndex[i][0]] - 1; j++)
			{
				gPepMassSeq[gSequenceNum][j] = gPepMassSeq[newSeqIndex[i][0]][j];
				testMass += gPepMassSeq[gSequenceNum][j];
			}
			massDiff = gPepMassSeq[newSeqIndex[i][1]][0] - testMass;
			gPepMassSeq[gSequenceNum][j] = massDiff;
			j++;
			for(k = 1; k < gPepLength[newSeqIndex[i][1]]; k++)
			{
				gPepMassSeq[gSequenceNum][j] = gPepMassSeq[newSeqIndex[i][1]][k];
				j++;
			}
			gPepLength[gSequenceNum] = j;
			gSequenceNum++;
		}
	}

	return;
}

/***************************FleshOutSequences********************************************************************
*
*	First make list of unsequenced masses (no repeats).  Then for each unsequenced mass start to find amino 
*	acid combinations that match, but force a K or R if its a C-terminal mass.  For each combination, 
*	make all possible sequences (leaving K or R at the C-term), and score them using a simple y/b score.
*	The best score wins and replaces the unsequenced mass.
*/

void FleshOutSequenceEnds(struct MSData *firstMassPtr)
{
	INT_4 i, j, k, l, mass;
	INT_4 cTermMassNum, nTermMassNum, *cTermMasses, *nTermMasses;
	INT_4 *sequenceToAdd, *sequenceToAppend;
	char  cTerm;
	
	
	/*Assign space to arrays*/
	cTermMasses 			= (int *) malloc(MAX_SEQUENCES * 2 * sizeof(INT_4));
	if(cTermMasses == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}
	
	nTermMasses 			= (int *) malloc(MAX_SEQUENCES * 2 * sizeof(INT_4));
	if(nTermMasses == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}
	
	sequenceToAdd 			= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(sequenceToAdd == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}
	
	sequenceToAppend 		= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(sequenceToAppend == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}
	
	/*Initialize*/
	cTermMassNum = 0;
	nTermMassNum = 0;
	for(i = 0; i < MAX_SEQUENCES*2; i++)
	{
		cTermMasses[i] = 0;
		nTermMasses[i] = 0;
	}
	
	/*Make a list of C-terminal unsequenced masses*/
	GetCTermMasses(&cTermMassNum, cTermMasses);
	
	/*Make a list of N-terminal unsequenced masses*/
	GetNTermMasses(&nTermMassNum, nTermMasses);
	
	/*Make a mass ordered list of monoisotopic amino acids */
	MakeAAArray();
	
	/*Find best sequence for each n-terminal mass*/
	for(i = 0; i < nTermMassNum; i++)
	{
		mass = nTermMasses[i];
		sequenceToAdd[0] = 0;	/*initialize*/
		GetBestNtermSeq(sequenceToAdd, mass, firstMassPtr);
		
		/*Now stick this bit of sequence onto the appropriate full sequences (replacement)*/
		if(sequenceToAdd[0] != 0)
		{
			for(j = 0; j < gSequenceNum; j++)
			{
				if(mass >= gPepMassSeq[j][0] - gParam.fragmentErr &&
					mass <= gPepMassSeq[j][0] + gParam.fragmentErr)
				{
					l = 0; 
					while(sequenceToAdd[l] != 0 && l < MAX_PEPTIDE_LENGTH)
					{
						sequenceToAppend[l] = sequenceToAdd[l];
						l++;
					}

					for(k = 1; k < gPepLength[j]; k++)	/*don't add k=0, since thats the unsequenced mass*/
					{
						if(l < MAX_PEPTIDE_LENGTH)
						{
							sequenceToAppend[l] = gPepMassSeq[j][k];
							l++;
						}
					}
					if(l < MAX_PEPTIDE_LENGTH)
					{
						gPepLength[j] = l;
						for(k = 0; k < gPepLength[j]; k++)
						{
							gPepMassSeq[j][k] = sequenceToAppend[k];
						}
					}
				}
			}
		}
	}
		
	/*Find best sequence for each c-terminal mass, assuming K or R at C-terminus*/
	
	/*Figure out if there is a C-terminal Lys, Arg, or Both*/
	cTerm = CheckCterm(firstMassPtr);
	
	/*Now proceed*/
	for(i = 0; i < cTermMassNum; i++)
	{
		mass = cTermMasses[i];
		sequenceToAdd[0] = 0;	/*initialize*/
		GetBestCtermSeq(sequenceToAdd, mass, cTerm, firstMassPtr);
		
		/*Now stick this bit of sequence onto the appropriate full sequences (replacement)*/
		if(sequenceToAdd[0] != 0)
		{
			for(j = 0; j < gSequenceNum; j++)
			{
				if(mass >= gPepMassSeq[j][gPepLength[j] - 1] - gParam.fragmentErr &&
					mass <= gPepMassSeq[j][gPepLength[j] - 1] + gParam.fragmentErr)
				{
					l = 0;
					while(sequenceToAdd[l] != 0)
					{
						l++;
					}
					if(l < MAX_PEPTIDE_LENGTH)
					{
						l = 0;
						k = gPepLength[j] - 1;
						while(sequenceToAdd[l] != 0)
						{
							gPepMassSeq[j][k] = sequenceToAdd[l];
							k++;
							l++;
						}
						gPepLength[j] = k;
					}
				}
			}
		}
	}
	
	/*free the arrays*/
	free(cTermMasses);
	free(nTermMasses);
	free(sequenceToAdd);
	free(sequenceToAppend);
	
	return;
}

/**********************************GetBestCtermSeq******************************************
*
*	For each input mass of unsequenced C-terminus, find all random sequences that fit the mass.
*	To limit the computation time, only masses less than an upper limit are examined.  Sequences
*	that fit the mass are scored according to how many y and b ions are matched.  A C-terminal
*	Arg and/or Lys is assumed for tryptic peptides.
*/

void	GetBestCtermSeq(INT_4 *sequenceToAdd, INT_4 mass, char cTerm, struct MSData *firstMassPtr)
{
	INT_4 	maxResidues, minResidues, residueNum, *sequence, i, j, k;	
	INT_4	ratchetMass, nominalMass, loopNumber, newMass;
	REAL_4	testNum, position, score, bestScore;
	REAL_4	massLimit = 747;	/*largest bit of unsequenced mass to be examined*/
	char	cTermAA;
	
	if((clock() - gParam.startTicks)/ CLOCKS_PER_SEC > 45)
	{
		massLimit = 600;
	}
	if(mass > massLimit * gMultiplier)
	{
		sequenceToAdd[0] = 0;	/*signal that nothing was found*/
		return;	/*only try for short bits of mass*/
	}
		
	/*Assign space to arrays*/
	sequence 			= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(sequence == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}

	/*Initialize*/
	for(i = 0; i < MAX_PEPTIDE_LENGTH; i++)
	{
		sequence[i] 		= 0;	/*this is the array used to find random sequences*/
		sequenceToAdd[i] 	= 0;	/*this is the final best sequence to add for this particular mass */
	}
	bestScore = 0;
	
	
	/*Create a loop that considers Arg, Lys, or both at C-terminus*/
	if(cTerm == 'B')
	{
		loopNumber = 2;	/*both y1 for Arg and Lys were found*/
	}
	else
	{
		loopNumber = 1;	/*only y1 for either Arg or Lys were found, *or* its not a tryptic peptide*/
	}
	
	for(k = 0; k < loopNumber; k++)
	{
		if(loopNumber == 2)
		{
			if(k == 0)
			{
				newMass = mass - gMonoMass_x100[K];
				cTermAA = 'K';	/*this denotes the particular C-term amino acid this time through the k loop*/
			}
			else
			{
				newMass = mass - gMonoMass_x100[R];
				cTermAA = 'R';
			}
		}
		else
		{
			if(cTerm == 'K')
			{
				newMass = mass - gMonoMass_x100[K];
				cTermAA = 'K';
			}
			else if(cTerm == 'R')
			{
				newMass = mass - gMonoMass_x100[R];
				cTermAA = 'R';
			}
			else
			{
				newMass = mass;
				cTermAA = 'N';
			}
		}
		

		/*Find max and min number of residues*/
		testNum 	= (REAL_4)newMass / gMonoMass_x100[G];
		maxResidues = testNum;
		if(maxResidues < 2)
		{
			return;	/*presumably a single amino acid would have been found already*/
		}
		if(maxResidues > MAX_PEPTIDE_LENGTH)
		{	
			maxResidues = MAX_PEPTIDE_LENGTH;
		}
		
		testNum 	= (REAL_4)newMass / gMonoMass_x100[W];
		testNum		= testNum + 1;
		minResidues	= testNum;
		if(minResidues < 2)
		{
			minResidues = 2;	/*presumably a single amino acid would have been found already*/
		}
		
		nominalMass = (REAL_4)newMass / gMultiplier + 0.5;	/*use nominal masses now*/
		
		/*start searching for sequences*/
		for(residueNum = minResidues; residueNum <= maxResidues; residueNum++)
		{
			/*Initialize each time the sequence length changes*/
			ratchetMass = 0;
			for(i = 0; i < residueNum; i++)
			{
				sequence[i] = 0;
				ratchetMass += gAAArray[0];
			}
			sequence[0] = -1;	/*The first time through Ratchet moves this to zero*/
			position	= 0;
			
			/*Ratchet produces a new sequence until all have been done, at which point it returns a NULL*/
			while(ratchetMass != 0)
			{
				/*Ratchets through all possible sequences, and returns mass of sequence*/
				ratchetMass = RatchetHaggis(sequence, residueNum, position, ratchetMass, nominalMass);
				
			
				if(ratchetMass == nominalMass)
				{
					/*Score the sequence*/
					if(cTermAA == 'K')
					{
						sequence[residueNum] = gCTermKIndex;
						score = CSequenceScore(mass, sequence, residueNum + 1, firstMassPtr);
					}
					else if(cTermAA == 'R')
					{
						sequence[residueNum] = gCTermRIndex;
						score = CSequenceScore(mass, sequence, residueNum + 1, firstMassPtr);
					}
					else
					{
						score = CSequenceScore(mass, sequence, residueNum, firstMassPtr);
					}
					
					/*Check if this is the highest scoring sequence so far (save it)*/
					if(score > bestScore)
					{
						bestScore = score;
						if(cTermAA == 'K' || cTermAA == 'R')
						{
							for(i = 0; i < residueNum + 1; i++)
							{
								sequenceToAdd[i] = gAAArray[sequence[i]];	/*puts nominal masses into array*/
							}
							for(i = residueNum + 1; i < MAX_PEPTIDE_LENGTH; i++)
							{
								sequenceToAdd[i] = 0;	/*backfill*/
							}
						}
						else
						{
							for(i = 0; i < residueNum; i++)
							{
								sequenceToAdd[i] = gAAArray[sequence[i]];	/*puts nominal masses into array*/
							}
							for(i = residueNum; i < MAX_PEPTIDE_LENGTH; i++)
							{
								sequenceToAdd[i] = 0;	/*backfill*/
							}
						}
					}
				}
			}
		}
	}
	
	/*Replace nominal masses in sequenceToAdd with monoisotopic masses*/
	i = 0;
	while(sequenceToAdd[i] != 0)
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(j != K && j != I)
			{
				if(sequenceToAdd[i] == gNomMass[j])
				{
					sequenceToAdd[i] = gMonoMass_x100[j];
					break;
				}
			}
		}
		i++;
	}

	/*free the arrays*/
	free(sequence);

	return;
}

/*********************************CSequenceScore*********************************************
*
*	Calculate a score for each candidate c-terminal sequence.
*/
REAL_4	CSequenceScore(INT_4 mass, INT_4 *sequence, INT_4 residueNum, struct MSData *firstMassPtr)
{
	REAL_4	score, precursor;
	INT_4	bIons[MAX_PEPTIDE_LENGTH], yIons[MAX_PEPTIDE_LENGTH];
	INT_4	i, j, maxCharge;
	struct MSData *currPtr;
	
	/*Initialize*/
	score = 0;
	if(gParam.chargeState > 2)
	{
		maxCharge = 2;
	}
	else
	{
		maxCharge = 1;
	}
	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;
	
	/*Calculate b ions at different charge states and search for them*/
	
	for(j = 1; j <= maxCharge; j++)
	{
	
		/*Calculate b ions*/
		bIons[0] = (gParam.peptideMW - gParam.modifiedCTerm - mass + (j-1) * gElementMass_x100[HYDROGEN]) / j;
		
		for(i = 1; i < residueNum; i++)
		{
			bIons[i] = bIons[i-1] + gAAMonoArray[sequence[i-1]] / j;
		}
		
		/*Look for the b ions*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= bIons[1] - gParam.fragmentErr)	/*skip over the b1 ions*/
			{
				if(currPtr->mOverZ > bIons[residueNum - 1] + gParam.fragmentErr)
				{
					break;	/*stop looking above the highest mass b ion*/
				}
				for(i = 1; i < residueNum; i++)	/*the highest b ion should not count*/
				{
					if(currPtr->mOverZ >= bIons[i] - gParam.fragmentErr &&
						currPtr->mOverZ <= bIons[i] + gParam.fragmentErr)
					{
						if(bIons[i] > 350 * gMultiplier *(j - 1))
						{
							if(bIons[i] < precursor || gParam.fragmentPattern == 'L')
							{
								score += currPtr->intensity;
							}
						}
					}
				}
			}
			currPtr = currPtr->next;
		}
	}
	
	/*Calculate y ions at different charge states and search for them*/
	for(j = 1; j <= maxCharge; j++)
	{
	
		/*Calculate y ions*/
		yIons[0] = (gParam.modifiedCTerm + 2 * gElementMass_x100[HYDROGEN] + mass 
					+ (j-1) * gElementMass_x100[HYDROGEN]) / j;	
		for(i = 1; i < residueNum; i++)
		{
			yIons[i] = yIons[i-1] - gAAMonoArray[sequence[i - 1]] / j;
		}
		
		/*Look for the y ions*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= yIons[residueNum - 1] - gParam.fragmentErr)	/*skip over lower masses*/
			{
				if(currPtr->mOverZ > yIons[1] + gParam.fragmentErr)
				{
					break;	/*stop looking above the highest mass y ion*/
				}
				for(i = 1; i <= residueNum; i++)	/*the highest mass y ion should not count*/
				{
					if(currPtr->mOverZ >= yIons[i] - gParam.fragmentErr &&
						currPtr->mOverZ <= yIons[i] + gParam.fragmentErr)
					{
						if(yIons[i] > 350 * gMultiplier *(j - 1))
						{
							score += currPtr->intensity;
						}
					}
				}
			}
			currPtr = currPtr->next;
		}
	}
	
	

	return(score);
}


/**********************************CheckCterm**********************************************
*
*
*/
char	CheckCterm(struct MSData *firstMassPtr)
{
	INT_4	yLys, yArg;
	char	cTerm;
	BOOLEAN	yLysFound, yArgFound;
	struct MSData *currPtr;
	
	/*Intialize*/
	yLys = gMonoMass_x100[K] + 2 * gElementMass_x100[HYDROGEN] + gParam.modifiedCTerm;
	yArg = gMonoMass_x100[R] + 2 * gElementMass_x100[HYDROGEN] + gParam.modifiedCTerm;
	yLysFound = FALSE;
	yArgFound = FALSE;
	
	/*Look for the y1 ions*/
	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ >= yLys - gParam.fragmentErr)
		{
			if(currPtr->mOverZ > yArg + gParam.fragmentErr)
			{
				break;
			}
			if(currPtr->mOverZ >= yLys - gParam.fragmentErr &&
				currPtr->mOverZ <= yLys + gParam.fragmentErr)
			{
				yLysFound = TRUE;
			}
			if(currPtr->mOverZ >= yArg - gParam.fragmentErr &&
				currPtr->mOverZ <= yArg + gParam.fragmentErr)
			{
				yArgFound = TRUE;
			}
		}
		currPtr = currPtr->next;
	}
	
	/*Now decide what to report back*/
	if(yArgFound && !yLysFound)
	{
		cTerm = 'R';
	}
	else if(yLysFound && !yArgFound)
	{
		cTerm = 'K';
	}
	else
	{
		cTerm = 'B';	/*neither or both were found, so its ambiguous*/
	}
	
	if(gParam.proteolysis != 'T')
	{
		cTerm = 'N';	/*its not a tryptic cleavage*/
	}

	return(cTerm);
}

/***********************************MakeAAArray********************************************
*
*
*/
void	MakeAAArray(void)
{
	INT_4 i, j, smallestNumber, smallestNumberIndex;
	BOOLEAN keep;
	
	/*Intialize*/
	for(i = 0; i < AMINO_ACID_NUMBER; i++)
	{
		gAAArray[i] = 0;
	}
	gAANum = 0;
	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(i != Q)
		{
			smallestNumber = 100000000;
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(j != Q)
				{
					if(gMonoMass[j] < smallestNumber && gMonoMass[j] > 0)
					{
						smallestNumberIndex = j;
						smallestNumber = gMonoMass[j];
					}
				}
			}
			gMonoMass[smallestNumberIndex] *= -1;
			keep = TRUE;
			for(j = 0; j < gAANum; j++)
			{
				if(smallestNumber == gAAArray[j])
				{
					keep = FALSE;
					break;
				}
			}
			if(keep)
			{
				gAAArray[gAANum] = smallestNumber;
				gAAMonoArray[gAANum] = gMonoMass_x100[smallestNumberIndex];
				if(smallestNumberIndex == K)
				{
					gCTermKIndex = gAANum;
				}
				if(smallestNumberIndex == R)
				{
					gCTermRIndex = gAANum;
				}
				gAANum++;
			}
		}
	}
	
	/*Set gMonoMass back to positive numbers*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(gMonoMass[i] < 0)
		{
			gMonoMass[i] *= -1;
		}
	}
	
	gMassRange = gAAArray[gAANum - 1] - gAAArray[0];
	
	return;
}

/**********************************GetBestNtermSeq******************************************
*
*	For each input mass of unsequenced N-terminus, find all random sequences that fit the mass.
*	To limit the computation time, only masses less than an upper limit are examined.  Sequences
*	that fit the mass are scored according to how many y and b ions are matched.
*/

void	GetBestNtermSeq(INT_4 *sequenceToAdd, INT_4 mass, struct MSData *firstMassPtr)
{
	INT_4 	maxResidues, minResidues, residueNum, *sequence, i, j;	
	INT_4	ratchetMass, nominalMass;
	REAL_4	testNum, position, score, bestScore;
	REAL_4	massLimit = 600;	/*largest bit of unsequenced mass to be examined*/
	
	
	if(mass > massLimit * gMultiplier)
	{
		sequenceToAdd[0] = 0;	/*signal that nothing was found*/
		return;	/*only try for short bits of mass*/
	}
	
	/*Assign space to arrays*/
	sequence 			= (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(sequence == NULL)
	{
		printf("Haggis:FleshOutSequences memory error");
		exit(1);
	}

	/*Initialize*/
	for(i = 0; i < MAX_PEPTIDE_LENGTH; i++)
	{
		sequence[i] 		= 0;	/*this is the array used to find random sequences*/
		sequenceToAdd[i] 	= 0;	/*this is the final best sequence to add for this particular mass */
	}
	bestScore = 0;

	/*Find max and min number of residues*/
	testNum 	= (REAL_4)mass / gMonoMass_x100[G];
	maxResidues = testNum;
	if(maxResidues < 2)
	{
		return;	/*presumably a single amino acid would have been found already*/
	}
	if(maxResidues > MAX_PEPTIDE_LENGTH)
	{
		maxResidues = MAX_PEPTIDE_LENGTH;	/*this should never happen...*/
	}
	
	testNum 	= (REAL_4)mass / gMonoMass_x100[W];
	testNum		= testNum + 1;
	minResidues	= testNum;
	if(minResidues < 2)
	{
		minResidues = 2;	/*presumably a single amino acid would have been found already*/
	}
	
	nominalMass = mass / gMultiplier;	/*use nominal masses now*/
	
	for(residueNum = minResidues; residueNum <= maxResidues; residueNum++)
	{
		/*Initialize each time the sequence length changes*/
		ratchetMass = 0;
		for(i = 0; i < residueNum; i++)
		{
			sequence[i] = 0;
			ratchetMass += gAAArray[0];
		}
		sequence[0] = -1;	/*The first time through Ratchet moves this to zero*/
		position	= 0;
		
		/*Ratchet produces a new sequence until all have been done, at which point it returns a NULL*/
		while(ratchetMass != 0)
		{
			/*Ratchets through all possible sequences, and returns mass of sequence*/
			ratchetMass = RatchetHaggis(sequence, residueNum, position, ratchetMass, nominalMass);
			
		
			if(ratchetMass == nominalMass)
			{
				/*Score the sequence*/
				score = NSequenceScore(mass, sequence, residueNum, firstMassPtr);
				
				/*Check if this is the highest scoring sequence so far (save it)*/
				if(score > bestScore)
				{
					bestScore = score;
					for(i = 0; i < residueNum; i++)
					{
						sequenceToAdd[i] = gAAArray[sequence[i]];	/*puts nominal masses into array*/
					}
					for(i = residueNum; i < MAX_PEPTIDE_LENGTH; i++)
					{
						sequenceToAdd[i] = 0;	/*backfill*/
					}
				}
			}
		}
	}
	
	/*Replace nominal masses in sequenceToAdd with monoisotopic masses*/
	i = 0;
	while(sequenceToAdd[i] != 0)
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(j != K && j != I)
			{
				if(sequenceToAdd[i] == gNomMass[j] && i < MAX_PEPTIDE_LENGTH)
				{
					sequenceToAdd[i] = gMonoMass_x100[j];
					break;
				}
			}
		}
		i++;
	}

	/*free the arrays*/
	free(sequence);

	return;
}

/*********************************NSequenceScore*********************************************
*
*
*/
REAL_4	NSequenceScore(INT_4 mass, INT_4 *sequence, INT_4 residueNum, struct MSData *firstMassPtr)
{
	REAL_4	score, precursor;
	INT_4	bIons[MAX_PEPTIDE_LENGTH], yIons[MAX_PEPTIDE_LENGTH];
	INT_4	i, j, maxCharge;
	struct MSData *currPtr;
	
	/*Initialize*/
	score = 0;
	if(gParam.chargeState > 2)
	{
		maxCharge = 2;
	}
	else
	{
		maxCharge = 1;
	}
	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass_x100[HYDROGEN]) / gParam.chargeState;
	
	/*Calculate b ions at different charge states and search for them*/
	for(j = 1; j <= maxCharge; j++)
	{
	
		/*Calculate b ions*/
		bIons[0] = (gParam.modifiedNTerm + gAAMonoArray[sequence[0]] + (j-1) * gElementMass_x100[HYDROGEN]) / j;
		
		for(i = 1; i < residueNum; i++)
		{
			bIons[i] = bIons[i-1] + gAAMonoArray[sequence[i]] / j;
		}
		
		/*Look for the b ions*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= bIons[1] - gParam.fragmentErr)	/*skip over the b1 ions*/
			{
				if(currPtr->mOverZ > bIons[residueNum - 1] + gParam.fragmentErr)
				{
					break;	/*stop looking above the highest mass b ion*/
				}
				for(i = 1; i < residueNum - 1; i++)	/*the highest b ion should not count*/
				{
					if(currPtr->mOverZ >= bIons[i] - gParam.fragmentErr &&
						currPtr->mOverZ <= bIons[i] + gParam.fragmentErr)
					{
						if(bIons[i] > 350 * gMultiplier *(j - 1))
						{
							if(bIons[i] < precursor || gParam.fragmentPattern == 'L')
							{
								score += currPtr->intensity;
							}
						}
					}
				}
			}
			currPtr = currPtr->next;
		}
	}
	
	/*Calculate y ions at different charge states and search for them*/
	for(j = 1; j <= maxCharge; j++)
	{
	
		/*Calculate y ions*/
		yIons[0] = (gParam.peptideMW - gParam.modifiedNTerm + 2 * gElementMass_x100[HYDROGEN]
					 + (j-1) * gElementMass_x100[HYDROGEN]) / j;
		
		for(i = 1; i <= residueNum; i++)
		{
			yIons[i] = yIons[i-1] - gAAMonoArray[sequence[i - 1]] / j;
		}
		
		/*Look for the y ions*/
		currPtr = firstMassPtr;
		while(currPtr != NULL)
		{
			if(currPtr->mOverZ >= yIons[residueNum] - gParam.fragmentErr)	/*skip over lower masses*/
			{
				if(currPtr->mOverZ > yIons[1] + gParam.fragmentErr)
				{
					break;	/*stop looking above the highest mass b ion*/
				}
				for(i = 1; i < residueNum; i++)	/*the lowest mass y ion should not count*/
				{
					if(currPtr->mOverZ >= yIons[i] - gParam.fragmentErr &&
						currPtr->mOverZ <= yIons[i] + gParam.fragmentErr)
					{
						if(yIons[i] > 350 * gMultiplier *(j - 1))
						{
							score += currPtr->intensity;
						}
					}
				}
			}
			currPtr = currPtr->next;
		}
	}

	return(score);
}

/*********************************RatchetHaggis**********************************************
*
*
*/

INT_4	RatchetHaggis(INT_4 *sequence, INT_4 residueNum, INT_4 position, INT_4 ratchetMass, INT_4 correctMass)
{
	INT_4 i;

	if(ratchetMass == 0)
	{
		return(0);	/*make sure that this recursive call unwinds itself*/
	}
	
	if(sequence[0] == -1)	/*just starting out*/
	{
		sequence[0] += 1;
	}
	else if(sequence[position] < gAANum - 1)	/*changing the right most amino acid*/
	{
		ratchetMass = ratchetMass - gAAArray[sequence[position]];
		sequence[position] += 1;
		ratchetMass += gAAArray[sequence[position]];
	}
	else	/*need to ratchet*/
	{
		if(sequence[position] + 1 >= gAANum && position < residueNum)
		{
			ratchetMass = ratchetMass - gAAArray[gAANum - 1];
			ratchetMass += gAAArray[0];
			sequence[position] = 0;
			position += 1;
			if(position >= residueNum)
			{
				return(0);
			}
			ratchetMass = RatchetHaggis(sequence, residueNum, position, ratchetMass, correctMass);
		}
		else
		{
			return(0);	/*you've reached the end of the road*/
		}
	}
	
	/*Try to skip if nowhere near the right mass*/
	if(/*sequence[0] == 0 &&*/ ratchetMass != 0 && position < residueNum)
	{
		if(ratchetMass + gMassRange < correctMass)	/*if the mass is nowhere near close enough*/
		{
			ratchetMass = ratchetMass - gAAArray[sequence[0]];
			ratchetMass += gAAArray[gAANum - 1];
			sequence[0] = gAANum - 1;
			position = 0;
			ratchetMass = RatchetHaggis(sequence, residueNum, position, ratchetMass, correctMass);
		}
		if(ratchetMass > correctMass)	/*if the mass is already to much*/
		{
			for(i = 0; i <= position; i++)
			{
				ratchetMass = ratchetMass - gAAArray[sequence[i]];
				ratchetMass += gAAArray[gAANum - 1];
				sequence[i] = gAANum - 1;
			}
			position = 0;
			ratchetMass = RatchetHaggis(sequence, residueNum, position, ratchetMass, correctMass);
		}
	}
	
	return(ratchetMass);
}


/**********************GetCTermMasses****************************************************
*
*	Make a list of c-terminal unsequenced masses.
*/

void	GetCTermMasses(INT_4 *cTermMassNum, INT_4 *cTermMasses)
{
	INT_4 	i, j;
	BOOLEAN	newMassTest;
	
	for(i = 0; i < gSequenceNum; i++)
	{
		newMassTest	= TRUE;
		for(j = 0; j < *cTermMassNum; j++)
		{
			if(gPepMassSeq[i][gPepLength[i] - 1] <= cTermMasses[j] + gParam.fragmentErr &&
				gPepMassSeq[i][gPepLength[i] - 1] >= cTermMasses[j] - gParam.fragmentErr)
			{
				newMassTest = FALSE;
				break;
			}
		}
		if(newMassTest && *cTermMassNum < MAX_SEQUENCES * 2)
		{
			cTermMasses[*cTermMassNum] = gPepMassSeq[i][gPepLength[i] - 1];
			*cTermMassNum += 1;
		}
	}
	
	return;
}

/**********************GetNTermMasses****************************************************
*
*	Make a list of n-terminal unsequenced masses.
*/

void	GetNTermMasses(INT_4 *nTermMassNum, INT_4 *nTermMasses)
{
	INT_4 	i, j;
	BOOLEAN	newMassTest;
	
	for(i = 0; i < gSequenceNum; i++)
	{
		newMassTest	= TRUE;
		for(j = 0; j < *nTermMassNum; j++)
		{
			if(gPepMassSeq[i][0] <= nTermMasses[j] + gParam.fragmentErr &&
				gPepMassSeq[i][0] >= nTermMasses[j] - gParam.fragmentErr)
			{
				newMassTest = FALSE;
				break;
			}
		}
		if(newMassTest && *nTermMassNum < MAX_SEQUENCES * 2)
		{
			nTermMasses[*nTermMassNum] = gPepMassSeq[i][0];
			*nTermMassNum += 1;
		}
	}
	
	return;
}

/******************LoadHaggisSequenceStruct********************************************
*
* 	LoadSequenceStruct puts the residue masses in the peptide[] field, peptide length,
*	score, nodeValue, gapNum, and nodeCorrection, in their fields.  Of
*	course, to do this the function finds some memory, and this value is returned as a pointer
*	to a struct of type Sequence (which now contains all of this data).
*/
struct Sequence *LoadHaggisSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
							INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection)
{
	struct Sequence *currPtr;
	INT_4 i;	

	currPtr = (struct Sequence *) malloc(sizeof(struct Sequence));
	if(currPtr == NULL)
	{
		printf("LoadSequenceStruct in Haggis:  Out of mammories");
		exit(1);
	}
	
	for(i = 0; i < peptideLength; i++)
	{
		currPtr->peptide[i] = peptide[i];
	}	
	currPtr->peptideLength = peptideLength;	
	currPtr->score = score;
	currPtr->gapNum = gapNum;
	currPtr->nodeValue = nodeValue;
	currPtr->nodeCorrection = nodeCorrection;
	currPtr->next = NULL;

	return(currPtr);
}


/****************LinkHaggisSubsequenceList**********************************************************
*
* 	This function adds a subsequence onto the existing linked list of structs of type Sequence.
*	It adds structs in order of their score fields, so that the first in the list has the
*	highest score and the last in the list has the lowest score.
*
*/

struct Sequence *LinkHaggisSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr)
{
	struct Sequence *lastPtr;
	char test = TRUE;
		
	if(firstPtr == NULL)	/*If this is the first struct of the list then do this.*/
		firstPtr = newPtr;
	else
	{
		/*Find the last sequence in the list*/
		lastPtr = firstPtr;
		while(lastPtr->next != NULL)
		{
			lastPtr = lastPtr->next;
		}
		lastPtr->next = newPtr;
	}
	
	return(firstPtr);
}

/***************************************GetSequenceOfResidues******************************
*
*	Convert sequence of nodes to sequence of amino acid residue masses.  If a mass does
*	not equal a known residue mass then save the mass as is.
*/
void	GetSequenceOfResidues(INT_4 *mass)
{
	INT_4 i, j, testMass, k, y1R, y1K, residueMass;
	BOOLEAN y1Found, weirdMass;
	
	/*Initialize*/
	for(i = 0; i < MAX_SEQUENCES*2; i++)
	{
		gMatchSeries[i] = -1;
	}

	//Create sequences
	for(i = 0; i < gSeqCount; i++)
	{
		//Assume y ions
		j = 0; 
		while(gSequenceNodes[i][j] != 0)
		{
			j++;	//find the end of the sequence
		}
		gPepLength[gSequenceNum] = j + 1;
		testMass = mass[gSequenceNodes[i][0]] - gParam.modifiedCTerm
														- gElementMass_x100[HYDROGEN]*2;
		gPepMassSeq[gSequenceNum][j] = ResidueMass(testMass);	/*provide mass from gMonoMass_x100 if possible*/
		j--;
		k = 1;
		while(j > 0)
		{
			testMass = mass[gSequenceNodes[i][k]] - mass[gSequenceNodes[i][k-1]];
			gPepMassSeq[gSequenceNum][j] = ResidueMass(testMass);	/*provide mass from gMonoMass_x100 if possible*/
			j--;
			k++;
		}
		testMass = gParam.peptideMW + gElementMass_x100[HYDROGEN] - mass[gSequenceNodes[i][k-1]];
		gPepMassSeq[gSequenceNum][0] = ResidueMass(testMass);	/*provide mass from gMonoMass_x100 if possible*/
		gSequenceNum++;
		
		//Assume b ions
		y1R = gMonoMass_x100[R] + gElementMass_x100[HYDROGEN]*2 + gParam.modifiedCTerm;
		y1K = gMonoMass_x100[K] + gElementMass_x100[HYDROGEN]*2 + gParam.modifiedCTerm;
		testMass = mass[gSequenceNodes[i][0]];
		if((testMass > y1K - gParam.fragmentErr &&
			testMass < y1K + gParam.fragmentErr) ||
			(testMass > y1R - gParam.fragmentErr &&
			testMass < y1R + gParam.fragmentErr))
		{
			y1Found = TRUE;
		}
		else
		{
			y1Found = FALSE;
		}
		
		if(!y1Found)	//Don't even store if a y1 ion for Arg or Lys found
		{
			testMass = mass[gSequenceNodes[i][0]] - gParam.modifiedNTerm;
			gPepMassSeq[gSequenceNum][0] = ResidueMass(testMass);
			gMatchSeries[gSequenceNum] = gSequenceNum - 1;	/*denotes that this is a b ion series repeat*/
			j = 1;
			while(gSequenceNodes[i][j] != 0)
			{
				testMass = mass[gSequenceNodes[i][j]] - mass[gSequenceNodes[i][j-1]];
				gPepMassSeq[gSequenceNum][j] = ResidueMass(testMass);
				j++;
			}
			testMass = gParam.peptideMW - gParam.modifiedCTerm - mass[gSequenceNodes[i][j-1]];
			gPepMassSeq[gSequenceNum][j] = ResidueMass(testMass);
			gPepLength[gSequenceNum] = j + 1;
			gSequenceNum++;
		}
	}
	
/*
*	Clean up the N- and C-terminal ends to make sure they do not contain odd masses, like "53", or something
*/
	for(i = 0; i < gSequenceNum; i++)
	{
		residueMass = gPepMassSeq[i][gPepLength[i] - 1];
	
		/*Check to see if the C-terminal mass makes any sense (think about adding this later)*/
		if(residueMass < 184.121 * gMultiplier - gParam.fragmentErr)
		{
			weirdMass = TRUE;
			if((residueMass < 174.064 * gMultiplier + gParam.fragmentErr &&
				residueMass > 174.064 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 172.048 * gMultiplier + gParam.fragmentErr &&
				residueMass > 172.048 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 172.085 * gMultiplier + gParam.fragmentErr &&
				residueMass > 172.085 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 171.064 * gMultiplier + gParam.fragmentErr &&
				residueMass > 171.064 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 170.106 * gMultiplier + gParam.fragmentErr &&
				residueMass > 170.106 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 168.090 * gMultiplier + gParam.fragmentErr &&
				residueMass > 168.090 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 158.069 * gMultiplier + gParam.fragmentErr &&
				residueMass > 158.069 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 156.09 * gMultiplier + gParam.fragmentErr &&
				residueMass > 156.09 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 154.074 * gMultiplier + gParam.fragmentErr &&
				residueMass > 154.074 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 144.055 * gMultiplier + gParam.fragmentErr &&
				residueMass > 144.055 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 142.074 * gMultiplier + gParam.fragmentErr &&
				residueMass > 142.074 * gMultiplier - gParam.fragmentErr))
			{
				weirdMass = FALSE;	/*these are low mass two aa residue masses*/
			}
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(residueMass < gMonoMass_x100[j] + gParam.fragmentErr &&
					residueMass > gMonoMass_x100[j] - gParam.fragmentErr)
				{
					weirdMass = FALSE;
					break;
				}
			}
		}
		else
		{
			weirdMass = FALSE;
		}
		if(weirdMass)	/*add the weird c-term mass to the penultimate c-term mass*/
		{
			gPepMassSeq[i][gPepLength[i] - 2] += gPepMassSeq[i][gPepLength[i] - 1];
			gPepLength[i] -= 1;
		}
	
		/*Check to see if the N-terminal mass makes any sense*/
		
		residueMass = gPepMassSeq[i][0];
		
		if(residueMass < 184.121 * gMultiplier - gParam.fragmentErr)
		{
			weirdMass = TRUE;
			if((residueMass < 174.064 * gMultiplier + gParam.fragmentErr &&
				residueMass > 174.064 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 172.048 * gMultiplier + gParam.fragmentErr &&
				residueMass > 172.048 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 172.085 * gMultiplier + gParam.fragmentErr &&
				residueMass > 172.085 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 171.064 * gMultiplier + gParam.fragmentErr &&
				residueMass > 171.064 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 170.106 * gMultiplier + gParam.fragmentErr &&
				residueMass > 170.106 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 168.090 * gMultiplier + gParam.fragmentErr &&
				residueMass > 168.090 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 158.069 * gMultiplier + gParam.fragmentErr &&
				residueMass > 158.069 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 156.09 * gMultiplier + gParam.fragmentErr &&
				residueMass > 156.09 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 154.074 * gMultiplier + gParam.fragmentErr &&
				residueMass > 154.074 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 144.055 * gMultiplier + gParam.fragmentErr &&
				residueMass > 144.055 * gMultiplier - gParam.fragmentErr) ||
				(residueMass < 142.074 * gMultiplier + gParam.fragmentErr &&
				residueMass > 142.074 * gMultiplier - gParam.fragmentErr))
			{
				weirdMass = FALSE;	/*these are low mass two aa residue masses*/
			}
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(residueMass < gMonoMass_x100[j] + gParam.fragmentErr &&
					residueMass > gMonoMass_x100[j] - gParam.fragmentErr)
				{
					weirdMass = FALSE;
					break;
				}
			}
		}
		else
		{
			weirdMass = FALSE;
		}
		
		if(weirdMass)	//add the weird c-term mass to the penultimate c-term mass
		{
			gPepMassSeq[i][0] += gPepMassSeq[i][1];
			for(j = 2; j < gPepLength[i]; j++)
			{
				gPepMassSeq[i][j - 1] = gPepMassSeq[i][j];
			}
			gPepLength[i] -= 1;
		}
	}

	return;
}



/******************************FindNodeSequences*********************************
*
*	Given the ways that the ions can be connected (gForwardNodeConnect and gBackwardNodeConnect),
*	find all possible pathways through the nodes.
*/
void	FindNodeSequences(INT_4 *mass)
{
	INT_4 i, j, k, inputIon;
	BOOLEAN test;
	
	/*Initialize*/
	gSeqCount = 0;
	for(i = 0; i < MAX_SEQUENCES; i++)
	{
		for(j = 0; j < MAX_ION_NUM; j++)
		{
			gSequenceNodes[i][j] = 0;
		}
	}

/*	Step through the nodes from low mass to high mass*/
	for(i = 1; i < gIonCount; i++)	/*This is the first node, which is incremented up to the top*/
	{
		if(gForwardNum[i] != 0 && gNotTooManySequences)	/*gNotTooManySequences signals that array is filled*/
		{
			gEdgeNum 		= 0;	/*Counts the edges as the tree is searched*/
			test 			= TRUE;	/*Becomes FALSE when all paths are searched from node i*/
			inputIon 		= i;
			/*Need to get forward/backward arrays with correct positive/negative values for starting over*/
			for(j = 0; j < gIonCount; j++)
			{
				for(k = 0; k < gBackwardNum[j]; k++)
				{
					if(gBackwardNodeConnect[j][k] > 0)
					{
						gBackwardNodeConnect[j][k] = -1 * gBackwardNodeConnect[j][k];
					}
				}
				for(k = 0; k < gForwardNum[j]; k++)
				{
					if(gForwardNodeConnect[j][k] < 0)
					{
						gForwardNodeConnect[j][k] = -1 * gForwardNodeConnect[j][k];
					}
				}
			
			}
			/*test is positive until the low mass terminal node has no remaining pathways*/
			while(test)
			{
				test = NodeStep(&inputIon, mass);
			}
		}
	}
	return;
}

/******************************SetupBackwardAndForwardNodes***************************
*
*/

void	SetupBackwardAndForwardNodes(INT_4 *mass)
{
	INT_4 i, j, k, l, massDiff;
	
	
	/*Initialize variables*/
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		gForwardNum[i] 	= 0;
		gBackwardNum[i] 	= 0;
		for(j = 0; j < AMINO_ACID_NUMBER; j++)
		{
			gForwardNodeConnect[i][j] 	= 0;
			gBackwardNodeConnect[i][j] 	= 0;
		}
	}
	
	/*	First assume fragment ions are all singly charged.*/

	for(i = 0; i < gIonCount; i++)
	{
		for(j = i; j < gIonCount; j++)
		{
			
			massDiff = mass[j] - mass[i];
			if(massDiff >= gMonoMass_x100[G] - gParam.fragmentErr &&
			massDiff <= gMonoMass_x100[W] + gParam.fragmentErr)
			{
				for(k = 0; k < gAminoAcidNumber; k++)
				{
					if(massDiff <= gMonoMass_x100[k] + gParam.fragmentErr &&
						massDiff>= gMonoMass_x100[k] - gParam.fragmentErr)
					{
						gBackwardNodeConnect[j][gBackwardNum[j]] = -i;
						gBackwardNum[j]++;
						gForwardNodeConnect[i][gForwardNum[i]] = j;
						gForwardNum[i]++;
						break;
					}
				}
			}
			
		}
	}
	
/*	Clean up the node connections.  If a connection is made that is comprised of a larger jump that could also be
	two smaller connections (ie, two Gly to Gly jumps versus a single Asn jump), then the larger connection is 
	eliminated.*/
	
	for(i = 0; i < gIonCount; i++)
	{
		if(gForwardNum[i] > 1)
		{
			for(j = 0; j < gForwardNum[i] - 1; j++)
			{
				for(l = j + 1; l < gForwardNum[i]; l++)
				{
					if(gForwardNodeConnect[i][l] > 0 && gForwardNodeConnect[i][j] > 0)
					{
						massDiff = mass[gForwardNodeConnect[i][l]] - mass[gForwardNodeConnect[i][j]];
						if(massDiff < gParam.fragmentErr)
						{
							gForwardNodeConnect[i][l] = 0;	/*G+V=R, for example*/
						}
						else
						{
							for(k = 0; k < gAminoAcidNumber; k++)
							{
								if(massDiff <= gMonoMass_x100[k] + gParam.fragmentErr &&
									massDiff>= gMonoMass_x100[k] - gParam.fragmentErr)
								{
									gForwardNodeConnect[i][l] = 0;
								}
							}
						}
					}
				}
			}
		}
	}
	
	/*Get rid of the zero value node forward connections.*/
	for(i = 0; i < gIonCount; i++)
	{
		for(j = 0; j < gForwardNum[i]; j++)
		{
			if(gForwardNodeConnect[i][j] == 0)
			{
				for(l = j; l < gForwardNum[i]; l++)
				{
					gForwardNodeConnect[i][l] = gForwardNodeConnect[i][l+1];
				}
				gForwardNum[i] -= 1;	
			}
		}
	}

	return;
}

/******************************LoadMassArrays***************************************
*
*	Fills mass array with ion masses.
*/

INT_4	*LoadMassArrays(INT_4 *mass, struct MSData *firstMassPtr, INT_4 charge)
{
	INT_4 i;
	struct MSData *currPtr;

	/*Initialize variables*/
	gIonCount		=	1;
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		mass[i] = 0;
	}
	
	/*Load mass array*/
	currPtr = firstMassPtr;
	if(currPtr == NULL)
	{
		printf("Problem in LoadMassArrays: Haggis");
		exit(1);	/*No ions!!*/
	}

	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > 400 * gMultiplier * (charge - 1))	/*make sure big enough to hold the charge*/
		{
			if(currPtr->mOverZ * charge < gParam.peptideMW - gMonoMass_x100[G])
			{
				if(currPtr->mOverZ >= gMonoMass_x100[K] + gParam.modifiedCTerm + 2*gElementMass_x100[HYDROGEN] 
						- gParam.fragmentErr)
				{
					mass[gIonCount] = currPtr->mOverZ * charge - ((charge-1)*gElementMass_x100[HYDROGEN]);
					gIonCount++;
				}
				if(gIonCount > MAX_ION_NUM)
				{
					printf("Problem in LoadMassArrays: Haggis");
					exit(1);	/*Too many ions; I'll exceed the array sizes.*/
				}
				}
		}
		currPtr = currPtr->next;
	}
	return(mass);
}
/*****************************ResidueMass***********************************************
*
*	Given an input INT_4 testMass, determine if it matches to a residue mass in the
*	array gMonoMass_x100, then the gMonoMass_x100 derived number is returned.  If its
*	not found, then the original is returned.
*
*/

INT_4	ResidueMass(INT_4 inputMass)
{
	INT_4	i, outputMass;
	BOOLEAN massFound = FALSE;
	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(inputMass <= gMonoMass_x100[i] + gParam.fragmentErr &&
			inputMass >= gMonoMass_x100[i] - gParam.fragmentErr)
		{
			massFound = TRUE;
			outputMass = gMonoMass_x100[i];
			break;
		}
	}
	if(!massFound)
	{
		outputMass = inputMass;
	}
	return(outputMass);
}
/*****************************NodeStep*************************************************
*
*	NodeStep is a recursive function that takes in a node position, and then moves forward
*	or backward in the graph, all the while counting the number of edges in the path and
*	keeping track of the longest path.  It returns either a TRUE or FALSE value.  The default
*	return is a TRUE value, until all paths have been followed from a given starting node.  If
*	all paths have been followed, then FALSE is returned to signal the end of the path-finding
*	from a particular starting node.  Given a node position, the function steps forward one edge if
*	the edge has not previously been followed, the gEdgeNum is incremented up one and compared to gMaxEdgeNum.
*	The edge forward is used is made impassable (given a negative value), but the backwards edge is made positive and
*	passable.  If a node position has no forward edges, then it follows a passable backward edge and the
*	that edge used to travel backwards is made impassable again.  The edge numbering is decremented and the
*	function calls itself again.  Eventually, the program gets back to the starting node, and no edges are passable 
*	from that node, and that is when the function returns a FALSE value to signal that its time to move on to
*	another starting node.
*/

BOOLEAN	NodeStep(INT_4 *nodeNum, INT_4 *nodeMass)
{
	
	BOOLEAN wayForward = FALSE;	/*assume that there are no edges towards high mass*/
	BOOLEAN wayBackward = FALSE;	/*assume that there is no edges leading backwards*/
	BOOLEAN keepGoing = TRUE;			/*becomes FALSE when no more paths to follow*/
	BOOLEAN failureTest;
	INT_4 i, j, newNode, oldNode;
	

	if(gIonCount > MAX_ION_NUM)	/*array boundary*/
	{
		printf("gIonCount > MAX_ION_NUM");
		exit(1);
	}
	if(gAminoAcidNumber > AMINO_ACID_NUMBER)
	{
		printf("gAminoAcidNumber > AMINO_ACID_NUMBER");
		exit(1);
	}
	
	/*8/20/03 to allow for using parts of the tree that were used before, but connected via a different node
	at the bottom; might need to get rid of this if it causes problems later*/
	
	for(i = *nodeNum + 1; i < gIonCount/*MAX_ION_NUM*/; i++)
	{
		for(j = 0; j < gForwardNum[*nodeNum]; j++)
		{
			if(gForwardNodeConnect[i][j] < 0)
			{
				gForwardNodeConnect[i][j] *= -1;
			}
		}
	}
	
	

/*	Figure out if I can go up in mass or not*/

	if(gForwardNum[*nodeNum] != 0 && *nodeNum < gIonCount)
	{
		for(i = 0; i < gForwardNum[*nodeNum]; i++)
		{
			if(gForwardNodeConnect[*nodeNum][i] > 0)
			{
				wayForward = TRUE;	/*There is an edge available that leads to higher mass nodes*/
				break;
			}
		}
	}
	
/* 	Now that its known if there is an edge to higher mass, or not, we can proceed.*/

	if(wayForward)
	{
		
		for(i = 0; i < gForwardNum[*nodeNum]; i++)
		{
			if(gForwardNodeConnect[*nodeNum][i] > 0 && i < gAminoAcidNumber && *nodeNum < gIonCount)
			{
				gEdgeNum++;
				newNode = gForwardNodeConnect[*nodeNum][i];
				/*if(gEdgeNum > gMaxEdgeNum)
				{
					gMaxEdgeNum = gEdgeNum;
					SaveLongestSequence(*nodeNum, newNode);
				}*/
				gForwardNodeConnect[*nodeNum][i] = -1 * gForwardNodeConnect[*nodeNum][i];
				failureTest = TRUE;
				for(j = 0; j < gBackwardNum[newNode]; j++)
				{
					if(*nodeNum == -1 * gBackwardNodeConnect[newNode][j])
					{
						gBackwardNodeConnect[newNode][j] = -1 * gBackwardNodeConnect[newNode][j];
						failureTest = FALSE;
						break;
					}
				}
				if(failureTest)
				{
					printf("Problem in function NodeStep");
					exit(1);
				}
				break;
			}
		}
		oldNode = *nodeNum;	/*debug*/
		*nodeNum = newNode;
		keepGoing = NodeStep(nodeNum, nodeMass);
		if(!keepGoing)
		{
			return(FALSE);
		}
	
	}
	else	/*need to back-track, if possible*/
	{
	
/*	Store the sequences that cannot be extended (in a global array)*/
		if(gSeqCount >= MAX_SEQUENCES)
		{
			gNotTooManySequences = FALSE;
			printf("Haggis had to quit early, because there were too many sequences found.\n");
			return(FALSE);	/*too many sequences, so stop now*/
		}
		else
		{
			StoreSeq(*nodeNum, nodeMass);	/*there is room for more sequences*/
		}
		
/*	Figure out if I can go down in mass or not*/

		if(gBackwardNum[*nodeNum] != 0)
		{
			for(i = 0; i <gBackwardNum[*nodeNum]; i++)
			{
				if(gBackwardNodeConnect[*nodeNum][i] > 0)
				{
					wayBackward = TRUE;	/*There is an edge available that leads to lower mass nodes*/
					break;
				}
			}
		}
		
		keepGoing = FALSE;	/*will keep going if there is an edge that leads backwards*/
		if(wayBackward)
		{
			for(i = 0; gBackwardNum[*nodeNum]; i++)
			{
				if(gBackwardNodeConnect[*nodeNum][i] > 0 && *nodeNum < gIonCount/*MAX_ION_NUM*/ && i < gAminoAcidNumber)
				{
					keepGoing = TRUE;
					gEdgeNum--;
					newNode = gBackwardNodeConnect[*nodeNum][i];
					gBackwardNodeConnect[*nodeNum][i] = -1 * gBackwardNodeConnect[*nodeNum][i];
					*nodeNum = newNode;
					break;
				}
			}
			keepGoing = NodeStep(nodeNum, nodeMass);
		}
		if(!keepGoing)
		{
			return(FALSE);
		}
		
	
	}
	return(keepGoing);
}


/**********************************StoreSeq************************************************************************
*
*	
*/
void	StoreSeq(INT_4 nodeNum, INT_4 *nodeMass)
{
	INT_4	i, j, k, index, gapNum, minEdgeNum;
	REAL_4 highMass, lowMass, testMass;
	BOOLEAN foundBottomNode = FALSE;
	BOOLEAN keepTheSeq = TRUE;
	BOOLEAN testForGap;
	
	/*Quit before running out of space*/
	if(gSeqCount >= MAX_SEQUENCES)
	{
		printf("Way too many sequences.");
		exit(1);
	}
	
	if(gEdgeNum >= MAX_ION_NUM)	/*check array boundaries*/
	{
		printf("Problem in StoreSeq");
		exit(1);
	}
	
	testMass = gParam.peptideMW / gMultiplier;	/*peptide mass*/
	testMass = testMass / AV_RESIDUE_MASS;	/*guess at the number of residues*/
	if(gLutefiskSequenceCount > 10000)
	{
		minEdgeNum = testMass / 2 + 0.5;
	}
	else if(gLutefiskSequenceCount > 1000)
	{
		minEdgeNum = testMass / 3 + 0.5;	/*need a series that covers a third of the sequence*/
	}
	else
	{
		minEdgeNum = testMass / 4 + 0.5;
	}
	if(minEdgeNum < 4)
		minEdgeNum = 4;	/*bottom limit*/

	if(gEdgeNum < minEdgeNum)	/*Anything with fewer than 4 edges is a useless sequence*/
	{
		return;
	}

/*	Initialize the gSequenceNodes*/
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		gSequenceNodes[gSeqCount][i] = 0;
	}

/*	Fill in the sequence*/
	gSequenceNodes[gSeqCount][gEdgeNum] = nodeNum;
	i = gEdgeNum;
	gapNum = 0;

	while(!foundBottomNode && i > 0)
	{
		foundBottomNode = TRUE;
		index = gSequenceNodes[gSeqCount][i];
		for(j = 0; j < gBackwardNum[index]; j++)
		{
			if(gBackwardNodeConnect[index][j] > 0)
			{
				foundBottomNode = FALSE;
				highMass = nodeMass[gSequenceNodes[gSeqCount][i]];
				i--;
				gSequenceNodes[gSeqCount][i] = gBackwardNodeConnect[index][j];
				lowMass = nodeMass[gSequenceNodes[gSeqCount][i]];
				testMass = highMass - lowMass;
				testForGap = TRUE;	/*start by assuming its a gap*/
				for(k = 0; k < gAminoAcidNumber; k++)
				{
					if(testMass < gMonoMass_x100[k] + gParam.fragmentErr &&
						testMass > gMonoMass_x100[k] - gParam.fragmentErr)
					{
						testForGap = FALSE;	/*its not a gap*/
						break;
					}
				}
				if(testForGap)
				{
					gapNum++;
				}
				break;
			}
		}
	}

/*	Test for too many gaps*/
	if(gapNum > gParam.maxGapNum)
	{
		keepTheSeq = FALSE;
	}
	else
	{
		keepTheSeq = TRUE;
	}

/*	Test to see if its a subset of a previous sequence*/
	if(keepTheSeq)
	{
		for(i = 0; i < gSeqCount; i++)
		{
			j = 0;
			while(gSequenceNodes[i][j] != 0 || j == 0)
			{
				if(gSequenceNodes[gSeqCount][0] == gSequenceNodes[i][j])
				{
					k = 0;
					while(gSequenceNodes[gSeqCount][k] != 0)
					{
						if(gSequenceNodes[gSeqCount][k] != gSequenceNodes[i][j+k])
						{
							break;	/*break out if they are not the same, then check below to see if it 
									reached the end*/
						}
						k++;
					}
					if(gSequenceNodes[gSeqCount][k] == 0)	/*if the end was reached, then its a subset*/
					{
						keepTheSeq = FALSE;
						break;
					}
				}
				j++;
			}
			if(!keepTheSeq)
			{
				break;
			}
		}
	}
	
	if(keepTheSeq)	/*if the sequence is kept, then the sequence counter is incremented up one*/
	{
		gSeqCount++;
	}
	else
	{
		i = 0;
		while(gSequenceNodes[gSeqCount][i] != 0)
		{
			gSequenceNodes[gSeqCount][i] = 0;	/*reinitialize to zero*/
			i++;
		}
	}
	return;
}