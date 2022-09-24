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
	6/96
	
	LutefiskSubseqMaker is a file containing the function SubsequenceMaker plus its associated
	functions.  It was written to be used as part of a program called "LutefiskXP", which
	is used to aid in the interpretation of CID data of peptides.  The general aim of this file
	(and the function SubsequenceMaker) is to use the graph sequenceNode to derive a list of 
	completed sequences that account for some of the CID data.  It uses a subsequencing 
	approach.  Rather than searching each limb and twig of the tree, I ignore those branches 
	that do not appear to lead to anything interesting.
*/

#include <stdio.h>
#include <stdlib.h>
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"

/*Globals for this file only.*/
struct Sequence *gFinalSequencePtr;
INT_4 gSubseqNum = 0;
INT_4 gAA1Max, gAA1Min, gAA2Max, gAA2Min, gAA1, gAA2;


char gCheckItOut = FALSE;	/*Equals TRUE if I want to follow the subsequence buildup,
							or FALSE if I want it to run in the normal mode.*/

INT_4 gCorrectSequence[50] = {	/*Used for checking if the correct sequence is remaining.*/
	20011, 8703, 11308, 9907, 5702, 25009, 9907, 17106, 15610, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0
};

/********************************** ClearLowestScore *******************************
*
*	This function finds the lowest scoring subsequence (the last one in the list),
*	and removes all subsequences of that score.
*/

void ClearLowestScore(struct Sequence *newSubsequencePtr, INT_4 *subseqNum, 
						INT_4 maxLastNode, INT_4 minLastNode)
{
	struct Sequence *currPtr, *previousPtr, *trashPtr, *lastPtr;
	INT_4 lowScore;
	
	/*Return if only zero or one new subsequence.*/
	if(newSubsequencePtr->next == NULL || newSubsequencePtr == NULL)
		return;
	
	/*Find the end of the list.*/
	currPtr = newSubsequencePtr->next;
	previousPtr = newSubsequencePtr;
	while(currPtr != NULL)
	{
		previousPtr = previousPtr->next;
		currPtr = currPtr->next;
	}
	
	/*Find the lowest score, which is from the subsequence at the end of the list.*/
	lowScore = previousPtr->score;
	
	/*If the first new subsequence also has the low score, then get rid of one subsequence.*/
	if(newSubsequencePtr->score == lowScore)
	{
		currPtr = newSubsequencePtr->next;
		previousPtr = newSubsequencePtr;
		if(currPtr->next == NULL)	/*Return if only three subsequences.*/
			return;
		while(currPtr->next != NULL)
		{
			previousPtr = previousPtr->next;
			currPtr = currPtr->next;
		}
		previousPtr->next = NULL;
		free(currPtr);	/*Free one subsequence so that the new one can be added on.*/
		*subseqNum = *subseqNum - 1;
	}
	/*If the first subsequence has a higher score than the low score, get rid of all
	subsequences with scores equal to the low score.*/
	else
	{
		currPtr = newSubsequencePtr->next;
		previousPtr = newSubsequencePtr;
		while(currPtr != NULL && currPtr->score != lowScore)
		{
			previousPtr = previousPtr->next;
			currPtr = currPtr->next;
		}
		previousPtr->next = NULL;	/*terminate the list just prior to the low score subseqs*/
		lastPtr = previousPtr;	/*remember the last good subsequence pointer*/
		
		trashPtr = currPtr;
	/*	if(gParam.proteolysis == 'T')	/*keep subsequences that are about to finish*/
	/*	{
			while(currPtr != NULL)
			{
				if(currPtr->nodeValue + gMonoMass_x100[R] >= minLastNode)
				{
					if(currPtr->nodeValue + gMonoMass_x100[R] <= maxLastNode)
					{
						currPtr->nodeValue = (currPtr->nodeValue) * -1;
					}
					if(currPtr->nodeValue + gMonoMass_x100[K] >= minLastNode)
					{
						if(currPtr->nodeValue + gMonoMass_x100[K] <= maxLastNode)
						{
							currPtr->nodeValue = (currPtr->nodeValue) * -1;
						}
					}
				}
				currPtr = currPtr->next;
			}
			currPtr = trashPtr;
			while(currPtr != NULL)
			{
				if(currPtr->nodeValue < 0)
				{
					currPtr->nodeValue = (currPtr->nodeValue) * -1;
					lastPtr->next = currPtr;
					lastPtr = currPtr;
					currPtr = currPtr->next;
				}
				else
				{
					trashPtr = currPtr;
					currPtr = currPtr->next;
					free(trashPtr);
				}
			}
		}		
		else
		{*/
			FreeSequenceStructs(trashPtr);	/*trash the remaining ones.*/
	/*	}
	*/	
		/*Count the subsequences*/
		*subseqNum = 0;
		currPtr = newSubsequencePtr;
		while(currPtr != NULL)
		{
			*subseqNum = *subseqNum + 1;
			currPtr = currPtr->next;
		}
	}

	return;
}

/**********************************LCQNterminalSubsequences****************************************
*
*	This function sets up the first batch of subsequences.  It starts with a single subsequence
*	containing the N-terminal group (usually hydrogen of mass 1) and then tries to connect with
*	nodes that are up to 559 units higher.  If any higher mass nodes can be connected to lower
*	mass nodes, then the higher mass nodes are eliminated because these will be incorporated
*	into the sequence as the subsequencing progresses.  A linked list of structs of type
*	Sequence is generated where the first struct in the list has the highest score, and the last
*	struct in the list has the lowest score. 
*/

struct Sequence *LCQNterminalSubsequences(SCHAR *sequenceNode, INT_4 maxLastNode, INT_4 lowSuperNode,
											INT_4 highSuperNode)
{
	struct Sequence *subsequencePtr;
	INT_4 i, j, k, m, testValue, nTerminus;
	INT_4 *extensions, *extScore, extNum;
	INT_4 *bestExtensions, *bestExtScore, bestExtNum;
	INT_4 highestExtensionScore;
	INT_4 threshold;
	INT_4 threeAALimit = gAminoAcidNumber*gAminoAcidNumber*gAminoAcidNumber;
	INT_4 score, *peptide, peptideLength, nodeValue, gapNum;
	INT_4 *threeAA, threeAANum, sum;
	INT_4 sameNum, averageExtension, *sameExtension;
	INT_2 nodeCorrection;
	char nTerminusPossible, duplicateFlag;
	char sameTest, doIt;
	
	sameExtension = (int *) malloc(gParam.fragmentErr * 20 * sizeof(INT_4));
	if(sameExtension == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}

	threeAA = (int *) malloc(threeAALimit * sizeof(INT_4));
	if(threeAA == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}
	
	extensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(extensions == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}
	extScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(extScore == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}
	bestExtensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(bestExtensions == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}
	bestExtScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(bestExtScore == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}


	peptide = (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4 ));
	if(peptide == NULL)
	{
		printf("LCQNterminalSubsequences:  Out of memory.");
		exit(1);
	}
	
/*	Fill in the masses for three amino acids.*/

	threeAANum = 0;
	for(i = 0; i < gAminoAcidNumber; i++)	/*Fill in the masses of the 3 AA extensions.*/
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			for(k = 0; k < gAminoAcidNumber; k++)
			{
				if(gGapList[i] != 0 && gGapList[j] != 0 && gGapList[k] != 0)
				{
					sum = gGapList[i] + gGapList[j] + gGapList[k];
					duplicateFlag = FALSE;
					for(m = 0; m < threeAANum; m++)
					{
						if(threeAA[m] == sum)
						{
							/* We already have this mass in threeAA so don't add it to the list. */
							duplicateFlag = TRUE;
							break;
						}
					}
					if(duplicateFlag == FALSE)
					{
						for(m = 0; m <= gGapListIndex; m++)
						{
							if(gGapList[m] == sum)
							{
							/* We already have this mass  so don't add it to the list. */
							duplicateFlag = TRUE;
							break;
							}
						}
					}
					if(!duplicateFlag && threeAANum < threeAALimit - 1)
					{
						threeAA[threeAANum] = sum;
						threeAANum++;
					}
				}
			}
		}
	}


/*	Now start finding the N-terminal pieces.*/

	extNum = 0;
	subsequencePtr = NULL;
			
	nTerminus = gParam.modifiedNTerm;
												
	
/*	Find the one, two, and three amino acid jumps from gGapList and threeAA.*/

	for(i = nTerminus + gMonoMass_x100[G]; i < gMonoMass_x100[W] * 3; i++)	/*step thru each node*/
	{
		if(sequenceNode[i] != 0 && i < gGraphLength)	/*ignore the nodes w/ zero evidence*/
		{
			nTerminusPossible = FALSE;	/*start assuming that this is not an extension*/
/*	Check for one and two amino acid extensions using the gGapList array.*/
			for(j = 0; j <= gGapListIndex; j++)	
			{
				if(gGapList[j] != 0)
				{
					if(i - nTerminus == gGapList[j])
					{
						nTerminusPossible = TRUE;	/*its a possible extension*/
						break;
					}
				}
			}
/*	Now check for three amino acid extensions if it wasn't a one or two aa extension.*/
			if(nTerminusPossible == FALSE)
			{
				for(j = 0; j < threeAANum; j++)
				{
					if(threeAA[j] != 0)
					{
					 	if(i - nTerminus == threeAA[j])
					 	{
					 		nTerminusPossible = TRUE;
					 		break;
					 	}
					}
				}
			}
			doIt = TRUE;	/*check for superNodes when sequencetag specified*/
			if(nTerminus < lowSuperNode && i > highSuperNode)
			{
				doIt = FALSE;
			}
			if(nTerminusPossible && i <= maxLastNode && doIt)	/*save as an extension?*/
			{
				extensions[extNum] = i - nTerminus;
				extScore[extNum] = sequenceNode[i];
				extNum++;
				if(extNum >= MAX_GAPLIST)
				{
					printf("LCQNTerminalSubsequences:  extNum >= MAX_GAPLIST\n");
					exit(1);
				}
			}
		}
	}
	
/*	Find extensions that are 1 node unit apart, and consolidate them.*/

	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] != 0)
		{
			sameNum = 0;
			averageExtension = extensions[i];
			for(j = 0; j < extNum; j++)
			{
				if(j != i && extScore[j] != 0)
				{
					sameTest = FALSE;
					for(k = 0; k < sameNum; k++)
					{
						if(extensions[sameExtension[k]] - extensions[j] == 1 ||
							extensions[j] - extensions[sameExtension[k]] == 1)
						{
							sameTest = TRUE;
						}
					}
					if(extensions[i] - extensions[j] == 1 ||
						extensions[j] - extensions[i] == 1 || sameTest)
					{
						sameExtension[sameNum] = j;
						averageExtension += extensions[j];
						sameNum++;
					}
				}
			}
			if(sameNum != 0)
			{
				averageExtension = ((float)averageExtension / (sameNum + 1)) + 0.5;	/*count the i extension
																			and round the value*/
				extensions[i] = averageExtension;	
				for(j = 0; j < sameNum; j++)
				{
					extScore[sameExtension[j]] = 0;
				}
			}
		}
	}			
	
/*	Get rid of the extensions that are within gParam.fragmentErr of each other.*/
	
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] != 0)
		{
			for(j = 0; j < extNum; j++)
			{
				if(extScore[j] != 0 && i != j)
				{
					if(extensions[i] <= extensions[j] + gParam.fragmentErr
						&& extensions[i] >= extensions[j] - gParam.fragmentErr)
					{
						if(extScore[i] >= extScore[j])
						{
							extScore[j] = 0;
						}
						else
						{
							extScore[i] = 0;
						}
					}
				}
			}
		}
	}
	
	
/*	If two extensions differ by the mass of an amino acid, then the higher mass one's 
	intensity is assigned a zero so that it is removed in the section below.*/
	
	if(extNum > 0)
	{
		for(i = 0; i < extNum; i++)
		{
			for(j = i + 1; j < extNum; j++)
			{
				testValue = extensions[j] - extensions[i];
				for(k = 0; k < gAminoAcidNumber; k++)
				{
					if(gGapList[k] != 0)
					{
						if(testValue <= gGapList[k] + gParam.fragmentErr &&
							testValue >= gGapList[k] - gParam.fragmentErr)
						{
							extScore[j] = 0;
						}
					}
				}
			}
		}
	}
	
	
/*
*	Now I need to find the best extensions, ie, the top maxExtNum of them and only if these
*	extensions are greater than the product of the highest score and extThresh.
*	 bestExtensions[MAX_GAPLIST], bestExtScore[MAX_GAPLIST], bestExtNum;
*/

	highestExtensionScore = extScore[0];	/*Find the highest extension score.*/
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] > highestExtensionScore)
		{
			highestExtensionScore = extScore[i];
		}
	}
	
	threshold = highestExtensionScore * gParam.extThresh; /*Set the extension score threshold.*/

/*
*	If a peptide tag has been entered and if the highestExtensionScore is over 100, then
*	that means that the highest scoring extension is a superNode.  In this case, don't
*	use an extension score threshold.
*/
	if(gParam.tagNMass != 0 && gParam.tagCMass != 0)
	{
		if(highestExtensionScore > 100)
		{
			threshold = 0;
		}
	}
	
	bestExtNum = 0;
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] >= threshold)
		{
			bestExtensions[bestExtNum] = extensions[i];
			bestExtScore[bestExtNum] = extScore[i];
			bestExtNum++;
			if(bestExtNum >= MAX_GAPLIST)
			{
				printf("LCQNTerminalSubsequences:  bestExtNum >= MAX_GAPLIST\n");
				exit(1);
			}
		}
	}
	
/*	Make sure there are enough subsequences allowed.*/
	if(bestExtNum > gParam.topSeqNum)
	{
		bestExtNum = gParam.topSeqNum;
	}
			
/*	Store this information in the linked list of Sequence structs.*/

	for(i = 0; i < bestExtNum; i++)
	{
		score = bestExtScore[i];
		peptide[0] = bestExtensions[i];
		peptideLength = 1;
		gapNum = 0;
		nodeValue = nTerminus + bestExtensions[i];
		j = gParam.modifiedNTerm * 10 + 0.5;
		k = gParam.modifiedNTerm + 0.5;
		k = k * 10;
		nodeCorrection = j - k;

		subsequencePtr = LinkSubsequenceList(subsequencePtr, LoadSequenceStruct(peptide,
																peptideLength,
																score, nodeValue, gapNum,
																nodeCorrection));
	}
	
/*	Free the arrays*/
	free(extensions);
	free(extScore);
	free(bestExtensions);
	free(bestExtScore);
	free(peptide);
	free(sameExtension);
	free(threeAA);
	
	return(subsequencePtr);
}


/*****************************amIHere******************************************************
*
*	This function is used to determine at each subsequence extension if the correct sequence
*	is present.  By setting gCheckItOut to be FALSE, this stuff is always skipped.  If 
*	gCheckItOut is TRUE, then this function is activated.
*
*/
void amIHere(INT_4 correctPeptideLength, struct Sequence *subsequencePtr)
{
	struct Sequence *currPtr, *correctPtr;
	INT_4 i, totalSubsequences, rank;
	INT_4 j = 0;
	INT_4 lowestScore;
	char test;
	
	if(subsequencePtr == NULL)
	{
		j++;	/*set debugger here*/
		j++;
		return;
	}

	
	/*gCheckItOut = FALSE;*/	/*If sequence is found then its made TRUE later.  This keeps it
							from running this function when the sequence drops out.*/
							
	currPtr = subsequencePtr;	/*Find the lowest score here.*/
	while(currPtr->next != NULL)
	{
		currPtr = currPtr->next;
	}
	lowestScore = currPtr->score;
	
	currPtr = subsequencePtr;
	while(currPtr != NULL)
	{
		if(currPtr->nodeValue == 13322)
		{
			j++;
			j++;
		}
		
		test = TRUE;
		for(i = 0; i < correctPeptideLength; i++)
		{
			if(currPtr->peptide[i] <= gCorrectSequence[i] - gParam.fragmentErr
				|| currPtr->peptide[i] >= gCorrectSequence[i] + gParam.fragmentErr)
			{
				test = FALSE;
				break;
			}
		}
		if(test)	/*If this is the correct subsequence, then..*/
		{
			gCheckItOut = TRUE;
			totalSubsequences = 0;
			rank = 1;
			correctPtr = currPtr;
			currPtr = subsequencePtr;
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
			j++;
			break;
		}
		currPtr = currPtr->next;
	}
	if(test != TRUE)
	{
		j++;	/*Stop here in the debugger for when it doesn't match anymore.*/
		j++;
	}
}

/*****************************CorrectMass**************************************************
*
*	This function figures out if the final sequence that is about to be stored, is in fact
*	of the correct mass.  It returns a FALSE if its not correct, and a TRUE if it is correct.
*/

char CorrectMass(INT_4 *peptide, INT_4 peptideLength, INT_4 *aaPresentMass)
{
	char correct = TRUE;
	char test, peptideChar[MAX_PEPTIDE_LENGTH];
	REAL_4 calcMass = 0;
	REAL_4 nMass;
	REAL_8 mToAFactor;
	INT_4 i, j, k, diff, tagSequenceMass;
	
	calcMass += gParam.modifiedNTerm;
	
	nMass = gParam.modifiedNTerm;
	
	calcMass += gParam.modifiedCTerm;
	
	for(i = 0; i < peptideLength; i++)	/*Figure out the nominal mass of the peptide.*/
	{
		calcMass += peptide[i];
	}
	
	
	if(calcMass >= gParam.monoToAv)	/*Convert from mono to average mass*/
	{
		mToAFactor = 0;
	}
	else
	{
		if(calcMass >= (gParam.monoToAv - gAvMonoTransition))
		{
			mToAFactor = (gParam.monoToAv - calcMass) / gAvMonoTransition;
		}
		else
		{
			mToAFactor = 1;
		}
	}
	mToAFactor = MONO_TO_AV - ((MONO_TO_AV - 1) * mToAFactor);
	
	if(calcMass > (gParam.monoToAv - gAvMonoTransition))		/*Convert from average to monoisotopic mass if necessary.*/
	{
		calcMass = calcMass * mToAFactor;
	}
	
/*	Now decide if the peptide mass is correct.*/
	if(calcMass <= (gParam.peptideMW + gParam.peptideErr) && calcMass >= (gParam.peptideMW - gParam.peptideErr))
	{
		correct = TRUE;
	}
	else
	{
		correct = FALSE;
	}
	
/*	Now decide if the sequence tag is present before storing it as a complete sequence.*/
	if(correct)
	{
		if(gParam.tagSequence[0] != '*')
		{
			test = TRUE;
			if(nMass >= (gParam.tagNMass - gParam.fragmentErr))	/*This is if the sequence tag starts at
													the N-terminus.*/
			{
				if(nMass <= (gParam.tagNMass + gParam.fragmentErr))
				{
					correct = TRUE;
					test = FALSE;
				}
			}

			if(test)
			{
				for(i = 0; i < peptideLength; i++)
				{
					nMass = nMass + (peptide[i]);
					if(nMass >= (gParam.tagNMass - gParam.fragmentErr))
					{
						if(nMass <= (gParam.tagNMass + gParam.fragmentErr))
						{
							correct = TRUE;
						}
						else
						{
							correct = FALSE;
						}
						break;
					}
				}
			}
		}
	}
/*	Now check to see if the amino acids that are supposed to be in the peptide are part of
	the sequence about to be strored away.*/
	if(correct)	/*If its not a correct sequence, then don't bother proving its even more
				incorrect.*/
	{	
		if(aaPresentMass[0] != -1)	/*Do this if there were any amino acids listed as being
									present in the sequence.*/
		{
/*	Identify 2 amino acid gaps in the sequence.*/
			for(i = 0; i < peptideLength; i++)
			{	
				test = TRUE;
				for(j = 0; j < gAminoAcidNumber; j++)
				{
					if(gGapList[j] != 0)
					{
						if(gGapList[j] == peptide[i])
						{
							test = FALSE;
							break;
						}
					}
				}
				if(test)
				{	
					peptideChar[i] = TRUE;	/*It is a 2 amino acid gap.*/
				}
				else
				{
					peptideChar[i] = FALSE;
				}
			}
/*	Now start comparing aaPresent and the sequence in 'peptide'.*/
			i = 0;
			while(aaPresentMass[i] != -1)
			{
				if(aaPresentMass[i] != 0)
				{
					test = TRUE;
					for(j = 0; j < peptideLength; j++)
					{
						if(peptide[j] == aaPresentMass[i])
						{
							test = FALSE;
							break;
						}
					}
					if(test)	/*I didn't find it as a 1 amino acid, but maybe its part of 2 aa.*/
					{
						for(j = 0; j < peptideLength; j++)
						{
							if(peptideChar[j])
							{
								diff = peptide[j] - aaPresentMass[i];
								for(k = 0; k < gAminoAcidNumber; k++)
								{
									if(gGapList[k] != 0)
									{
										if(diff >= gGapList[k] - gParam.fragmentErr
											&& diff <= gGapList[k] + gParam.fragmentErr)
										{
											test = FALSE;
											break;
										}
									}
								}
							}
						}
					}
					if(test)	/*I didn't find it in the sequence. Maybe its in the sequence tag.*/
					{
						j = 0;
						while(gParam.tagSequence[j] != 0)
						{
							tagSequenceMass = 0;
							for(k = 0; k < gAminoAcidNumber; k++)
							{
								if(gSingAA[k] == gParam.tagSequence[j])
								{
									if(gGapList[k] != 0)
									{
										tagSequenceMass = gGapList[k];
									}
									break;
								}
							}
							if(tagSequenceMass == aaPresentMass[i])
							{
								test = FALSE;
								break;
							}
							j++;
						}
					}
					if(test)
					{
						correct = FALSE;
						break;
					}
				}
			i++;
			}
		}
	}									

		
	return(correct);
}

/****************************AlterSubsequenceList**********************************************
*
* 	This function adds a subsequence onto the existing linked list of structs of type Sequence.
*	It adds structs in order of their score fields, so that the first in the list has the
*	highest score and the last in the list has the lowest score.
*
*/

struct Sequence *AlterSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr)
{
	struct Sequence *currPtr, *previousPtr;
	char test = TRUE;
	INT_4 i = 0;
	
	
	previousPtr = firstPtr;	/*This is to make sure currPtr and previousPtr are initialized, 
						so that in the end I can find the last element in the list to free.*/
						
						
	if(firstPtr->next != NULL)	/*Take care of situations where there is only one subsequence.*/
	{
		currPtr = firstPtr->next;
	}
	else
	{
		if(newPtr->score > firstPtr->score)
		{
			free(firstPtr);
			firstPtr = newPtr;
			return(firstPtr);
		}
		else
		{
			free(newPtr);
			return(firstPtr);
		}
	}
	
/*Much of this is the same as from LinkSubsequenceList.*/

	if(firstPtr == NULL) 	/*If this is the first struct of the list then do this.*/
	{
		firstPtr = newPtr;
		return(firstPtr);
	}	
	else
	{
		if(newPtr->score > firstPtr->score)	/*If the struct to be added has the best score then*/
		{
			newPtr->next = firstPtr;
			firstPtr = newPtr;
		}
		else	/*Otherwise, go find the position for the new struct (based on its score field.*/
		{
			while(currPtr->next != NULL)
			{
				if(newPtr->score > currPtr->score)	/*I found the place.*/
				{
					previousPtr->next = newPtr;
					newPtr->next = currPtr;
					test = FALSE;
					break;
				}
				previousPtr = currPtr;
				currPtr = currPtr->next;
			}
			if(test)	/*Ok, I didn't find the place, and I'm at the end of the linked list,
						so this must have the lowest score of all.  I'll leave the list alone
						and return.	This is different from LinkSubsequenceList, which added
						the new struct on to the end of the list.*/
			{
				free(newPtr);
				return(firstPtr);
			}
		}
	}
	
/*	Here's what's different from LinkSubsequenceList - I find the lowest score (or last in the
*	list) and I eliminate it.  I set the previousPtr's next field to NULL to make it the new
*	last element in the list, and I free the currPtr.
*/
	while(currPtr->next != NULL)
	{
		previousPtr = currPtr;
		currPtr = currPtr->next;
	}
	previousPtr->next = NULL;
	/*if(currPtr->nodeValue == 13322 && currPtr->score == 236)
	{
		i++;
		i++;
	} for debugging*/
		
	free(currPtr);
	
	return(firstPtr);
}
/**********************FreeSequenceStructs****************************************
*
* Used for freeing memory in a linked list.  Bob DuBose tells me its best to free 
* space in the reverse order
* that the space was malloc'ed.  This routine does that very thing.
*
*/

void FreeSequenceStructs(struct Sequence *s)
{
	struct Sequence *currPtr, *nextPtr;

     currPtr = s;
     while(currPtr != NULL)
     {
     	nextPtr = currPtr->next;
     	free(currPtr);
     	currPtr = nextPtr;
     }

    return;                             /* now unwind the recursion*/
}

/************************************* StoreSubsequences ********************************************
*
*  Store this information in the linked list of Sequence structs.  The values placed in the
*   extensionList are used to determine
*	values for the variables peptide[], score, peptideLength, gapNum, and nodeValue.  These
*	variables are passed to a few functions that are used to set up the linked list of structs
*	of type Sequence, which contain the next set of subsequences (newSubsequencePtr	
*/
struct Sequence *StoreSubsequences(
	struct Sequence *newSubsequencePtr,
	struct extension *extensionList,
	struct Sequence *currentSubsequence,
	INT_4 *lastNode, 
	INT_4 lastNodeNum,
	INT_4 maxLastNode, 
	INT_4 minLastNode,
	INT_4 *aaPresentMass, 
	INT_4 *seqNum,
	INT_4 *subseqNum, 
	SCHAR *sequenceNode)

{

	INT_4	i, j;				/* Loop index */
	BOOLEAN test;
	INT_4 *peptide, score, peptideLength, nodeValue, gapNum;	
	INT_2 nodeCorrection;
	
	
	peptide = (int *) calloc(MAX_PEPTIDE_LENGTH, sizeof(INT_4 ));
	if(peptide == NULL)
	{
		printf("StoreSubsequences:  Out of memory");
		exit(1);
	}


	/*Set up the peptide field so that it includes the previous values contained
	  in the peptide field of currPtr (the subsequence currently under investigation.*/
	for(i = 0; i < currentSubsequence->peptideLength; i++)	
	{
		peptide[i] = currentSubsequence->peptide[i];
	}
	
			
	i = 0;
	while(extensionList[i].mass > 0 &&  i < MAX_GAPLIST) 
	{						
		test = TRUE;	/*Becomes FALSE if the data is stored as a final sequence.*/
		score = currentSubsequence->score + extensionList[i].score;
		peptideLength = currentSubsequence->peptideLength + 1;
		if(peptideLength > MAX_PEPTIDE_LENGTH)
		{
			printf("StoreSubsequences: peptideLength > MAX_PEPTIDE_LENGTH\n");
			exit(1);
		}
		peptide[peptideLength - 1] = extensionList[i].mass;
		gapNum    = currentSubsequence->gapNum + extensionList[i].gapSize;
		nodeValue = currentSubsequence->nodeValue + extensionList[i].mass;
		nodeCorrection = currentSubsequence->nodeCorrection 
							+ extensionList[i].nodeCorrection;
		if(nodeCorrection >= 10)
		{
			nodeCorrection = nodeCorrection - 10;
			nodeValue = nodeValue + 1;
		}
		else if(nodeCorrection <= -10)
		{
			nodeCorrection = nodeCorrection + 10;
			nodeValue = nodeValue - 1;
		}
		
		if(nodeValue >= gAA1Min)	/*can be terminated by known C-terminal aa*/
		{
			if(nodeValue <= gAA1Max)
			{
				score = score + sequenceNode[maxLastNode];
				peptideLength++;
				if(peptideLength > MAX_PEPTIDE_LENGTH)
				{
					printf("StoreSubsequences: peptideLength > MAX_PEPTIDE_LENGTH\n");
					exit(1);
				}
				peptide[peptideLength - 1] = gAA1;
				nodeValue = nodeValue + gAA1;
			}
		}
		if(nodeValue >= gAA2Min)	/*can be terminated by known C-terminal aa*/
		{
			if(nodeValue <= gAA2Max)
			{
				score = score + sequenceNode[maxLastNode];
				peptideLength++;
				if(peptideLength > MAX_PEPTIDE_LENGTH)
				{
					printf("StoreSubsequences: peptideLength > MAX_PEPTIDE_LENGTH\n");
					exit(1);
				}
				peptide[peptideLength - 1] = gAA2;
				nodeValue = nodeValue + gAA2;
			}
		}
				
	
		/*Here's where the completed sequences are stored.*/
		if(nodeValue >= minLastNode && nodeValue <= maxLastNode)
		{
			for(j = 0; j < lastNodeNum; j++)
			{	
				if(nodeValue == lastNode[j])
				{
					test = FALSE;	/*Even if its not stored, I won't continue seq'ing this.*/
					
					if((gapNum <= gParam.maxGapNum) &&
							CorrectMass(peptide, peptideLength, aaPresentMass))
					{
						if(*seqNum < gParam.finalSeqNum)
						{
							gFinalSequencePtr = 
								LinkSubsequenceList(gFinalSequencePtr,
													LoadFinalSequenceStruct(peptide,
																	peptideLength,
																	score, nodeValue, 
																	gapNum, nodeCorrection));
							*seqNum = *seqNum + 1;
						}
						else
						{
							gFinalSequencePtr =
								AlterSubsequenceList(gFinalSequencePtr,
													LoadFinalSequenceStruct(peptide,
																	peptideLength,
																	score, nodeValue, 
																	gapNum, nodeCorrection));
						}
					}
				}
			}
		}			

		if(test && nodeValue < minLastNode)
		{
			if(gapNum <= gParam.maxGapNum)	/*Don't store if there are too many gaps.*/
			{
				if(*subseqNum < gParam.topSeqNum)	/*If there are not too many subsequences 
													  stored, then do this.*/
				{
					newSubsequencePtr = 
						LinkSubsequenceList(newSubsequencePtr, 
											LoadSequenceStruct(peptide, 
															peptideLength, 
															score, nodeValue, gapNum,
															nodeCorrection));
					*subseqNum = *subseqNum + 1;
				}
				else	/*If I have the max allowed number of subsequence, 
							then do this.*/
				{
					/*This is the new way, which removes all subsequences with the
					same low score*/
					/*ClearLowestScore(newSubsequencePtr, subseqNum, maxLastNode, 
									minLastNode);
					
					newSubsequencePtr = 
						LinkSubsequenceList(newSubsequencePtr, 
											LoadSequenceStruct(peptide, 
															peptideLength, 
															score, nodeValue, gapNum,
															nodeCorrection));
					*subseqNum = *subseqNum + 1;*/
					
					/*This is the old way which was to replace subsequences one at a time*/
					newSubsequencePtr = 
						AlterSubsequenceList(newSubsequencePtr, 
											LoadSequenceStruct(peptide, 
															peptideLength,
															score, nodeValue, gapNum,
															nodeCorrection));
				}

			}
		}
			
		i++;
	}
	free(peptide);
	
	return newSubsequencePtr;
}		

/************************************* ExtensionsSortDescend ********************************************
*  
*  qsort sorting subroutine used in SortExtensions().	
*/
int ExtensionsSortDescend(const void *n1, const void *n2) 
{
	struct extension *n3 = (struct extension *)n1;
	struct extension *n4 = (struct extension *)n2;

	if (n3->singleAAFLAG > n4->singleAAFLAG) {
		return -1;
	}	
	else if (n3->singleAAFLAG < n4->singleAAFLAG) {
		return 1;
	}
	else {
		if (n3->score > n4->score) {
			return -1;
		}
		else if (n3->score < n4->score) {
			return 1;
		}
		else return 0;
	}		
}


/************************************* SortExtensions ********************************************
*	
*/
struct extension *SortExtensions(struct extension *inExtensionList)

{
	INT_4 i;
	INT_4 extensionCount;
	INT_4 highestExtensionScore, threshold;
	INT_4 outIndex, lastOneInScore;
	struct extension *outExtensionList;

	outExtensionList = (extension *) calloc(MAX_GAPLIST, sizeof(struct extension)); /* Note: memory is zeroed */
	if(outExtensionList == NULL)
	{
		printf("SortExtensions:  Out of memory");
		exit(1);
	}
	
	/*Find the highest extension score.*/
	highestExtensionScore = inExtensionList[0].score;
	extensionCount = 0;
	i = 0;
	while (	inExtensionList[i].mass > 0 && i < MAX_GAPLIST ) 
	{
		if(inExtensionList[i].score > highestExtensionScore)
		{
			highestExtensionScore = inExtensionList[i].score;
		}
		i++;
		extensionCount++;
	}
	
	/*Set the extension score threshold.*/
	threshold = highestExtensionScore * gParam.extThresh; 	
/*	
*	If a peptide tag has been entered and if the highestExtensionScore is over 100, then 
*	that means that the highest scoring extension is to a superNode.  In this case, don't
*	use an extension score threshold.
*/
	if(gParam.tagNMass != 0 && gParam.tagCMass != 0)
	{
		if(highestExtensionScore > 100)
		{
			threshold = 0;
		}
	}	

	/*If the number of extensions found are less than the specified
	  maximum number of extensions per subsequence, then just  
	  put them in the outgoing list.*/
	if(extensionCount <= gParam.maxExtNum)	
	{
		outIndex = 0;
		for(i = 0; i < extensionCount; i++)
		{
			if(inExtensionList[i].score >= threshold && i < gParam.topSeqNum)
			{
				outExtensionList[outIndex] = inExtensionList[i];
				outIndex++;
				if(outIndex >= MAX_GAPLIST) {
					printf("FATAL ERROR in SortExtensions(): Overflowed the outExtensionList ?????");
					exit(1);
				}	
			}
		}
	}
	
	if(extensionCount > gParam.maxExtNum)	
	{
		/* There are too many extensions. Throw out the worst. */
		
		/* First, sort the extensions by whether they are singleAA or not and then by score. */
		qsort(inExtensionList, extensionCount, sizeof(struct extension), ExtensionsSortDescend);
		
		outIndex = 0;
		i = 0;
		while (i < extensionCount && outIndex < gParam.maxExtNum && outIndex < gParam.topSeqNum) {
			if (inExtensionList[i].score >= threshold)
			{
				outExtensionList[outIndex] = inExtensionList[i];
				outIndex++;
				if(outIndex >= MAX_GAPLIST) {
					printf("FATAL ERROR in SortExtensions(): Overflowed the outExtensionList ?????");
					exit(1);
				}	
			}
		
			i++;				
		}
		lastOneInScore = outExtensionList[outIndex - 1].score;
		/*now stick in the extensions with the same score as the lowest scoring extension*/
		while(i < extensionCount && outIndex < gParam.topSeqNum)
		{
			if(inExtensionList[i].score == lastOneInScore)
			{
				outExtensionList[outIndex] = inExtensionList[i];
				outIndex++;
				if(outIndex >= MAX_GAPLIST) {
					printf("FATAL ERROR in SortExtensions(): Overflowed the outExtensionList ?????");
					exit(1);
				}
			}
			i++;
		}
	}
	
	free(inExtensionList);
	
	return outExtensionList;
}
/************************************* Score2aaExtension ********************************************
*
*  Examine some special cases and then calculate the score for the 2 AA extension	
*/
struct extension Score2aaExtension(
	struct extension io2aaExtension,
	INT_4 startingNodeMass,
	INT_4 *oneEdgeNodes, 
	INT_4 oneEdgeNodesIndex,
	BOOLEAN oneEdgeNNode,
	char	sequenceNodeValue)
{
	INT_4 	i, j;   							/* Loop index */
	INT_4	massDiff;
	INT_4	endingNodeMass  = startingNodeMass + io2aaExtension.mass;
	BOOLEAN prolinePossible = FALSE;
/*	BOOLEAN glycinePossible = FALSE;*/
	BOOLEAN oneEdgeCNode    = FALSE;
	BOOLEAN precursorRegion = FALSE;
	REAL_4 precursor;
	INT_4 precursorMinErr, precursorPlusErr;

	precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN]) / gParam.chargeState;
	precursorPlusErr = precursor + 2 * gParam.fragmentErr + 0.5;
	precursorMinErr  = precursor - 2 * gParam.fragmentErr;


	massDiff = io2aaExtension.mass - gGapList[P];	/*Subtract the mass of proline.*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(gGapList[i] != 0)
		{
			if(massDiff <= gGapList[i] + gParam.fragmentErr &&
				massDiff >= gGapList[i] - gParam.fragmentErr)
			{
				prolinePossible = TRUE;
				break;
			}
		}
	}
	
	/*
	massDiff = io2aaExtension.mass - gGapList[G];	
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(gGapList[i] != 0)
		{
			if(massDiff == gGapList[i])
			{
				prolinePossible = true;
				break;
			}
		}
	}
	*/
	
	/*Find out if this node is a C-terminal one-edger.*/
	for(i = 0; i < oneEdgeNodesIndex; i++)
	{
		if(endingNodeMass == oneEdgeNodes[i])
		{
			oneEdgeCNode = TRUE;
			break;
		}
	}
	
	if(gParam.chargeState == 2)	/*For +2 precursors, if the 2 aa extension
								extends over the precursor ion, then don't
								attenuate the score so much.*/
	{
		if(startingNodeMass + gMonoMass_x100[G] <= precursor
		   && endingNodeMass - gMonoMass_x100[G] >= precursor)
		{
			for(i = 0; i < gAminoAcidNumber; i++)
			{
				if((startingNodeMass + gGapList[i] >= precursorMinErr)
				   && (startingNodeMass + gGapList[i] <= precursorPlusErr))
				{
					if(endingNodeMass - precursorMinErr >= gMonoMass_x100[G]
					   && endingNodeMass - precursorPlusErr <= gMonoMass_x100[W])
					{
						for(j = 0; j < gAminoAcidNumber; j++)
						{
							if(endingNodeMass - gGapList[j] >= precursorMinErr
							   && endingNodeMass - gGapList[j] <= precursorPlusErr)
							{
								precursorRegion = TRUE;
								break;
							}
						}
					}
					break;
				}
			}
		}
	}
	
/*	Now score the two amino acid extension depending on if its connecting two one-edge
*	nodes, one regular node and a one-edge node, or two regular nodes.
*/
	if(oneEdgeCNode && oneEdgeNNode)
	{
		/* Don't penalize as heavily for making a 2 aa jump that ends at the C-terminus. */
		io2aaExtension.score   = sequenceNodeValue * EDGE_EDGE_PENALTY;
		io2aaExtension.gapSize = 0;
	}
	else
	{
		if(prolinePossible)
		{
			io2aaExtension.score = sequenceNodeValue * PROLINE_PENALTY;
			io2aaExtension.gapSize = 0;
		}
		else
		{
			if(precursorRegion)
			{
				io2aaExtension.score = sequenceNodeValue * PRECURSOR_PENALTY;
			}
			else
			{
				/*if(glycinePossible)
				{
					io2aaExtension.score = sequenceNodeValue * GLYCINE_PENALTY;
				}
				else
				{*/
				if(oneEdgeCNode || oneEdgeNNode)
				{
					io2aaExtension.score = sequenceNodeValue * NODE_EDGE_PENALTY;
				}
				else
				{
					if(oneEdgeCNode == FALSE && oneEdgeNNode == FALSE)
					{
						io2aaExtension.score = sequenceNodeValue * NODE_NODE_PENALTY;
					}
				}	
				/*}*/
			}
		}
	}
	
	return io2aaExtension;
}		
/*************************************AddExtensions********************************************
*
*	This function adds single amino acid extensions onto the subsequences found in the linked
*	list starting w/ subsequencePtr.  As a result, a new linked list starting with 
*	newSubsequencePtr is generated.  When all of the extensions have been made and the list
*	starting with newSubsequencePtr is complete, then the list starting with subsequencePtr
*	is free'ed, and newSubsequencePtr is returned.  In addition to making one amino acid
*	jumps between nodes, this function performs two amino acid jumps for subsequences that 
*	cannot be extended any further, which leads to no penalty in the extension score.  If
*	a jump is made between a one-edge node and a regular node, then this extension score is
*	penalized by 0.75.  If a jump is made between two nodes that can be extended by single
*	amino acid jumps, then that extension score is penalized by 0.5.
*	The scores for the subsequences are derived from the array
*	sequenceNode, and are a summation of the scores for the nodes within a subsequence.  The
*	array gGapList and gGapListIndex describe the possible extensions that are allowed.
*	The array aaPresent is a character array listing in single letter code, the amino acids
*	known to be present.  If a completed sequence lacks one or more of these amino acids
*	then it is not stored.  extThresh is the extension threshold as described earlier, and
*	maxExtNum is the upper limit on the number of extensions allowed per subsequence.  The
*	maxGapNum is the upper limit on the number of two amino acid extensions allowed per
*	sequence (excluding the N-termninal two amino acids and certain C-terminal cases).  The
*	finalSeqNum is the upper limit on the number of completed sequences that will be stored
*	in the linked list of Sequence structs starting w/ gFinalSequencePtr.  The topSeqNum
*	is the maximum number of subsequences that are in the linked list of sequences (both
*	the list starting with subsequencePtr and newSubsequencePtr).
*/
struct Sequence *AddExtensions(struct Sequence *subsequencePtr, SCHAR *sequenceNode, 
						INT_4 *oneEdgeNodes, INT_4 oneEdgeNodesIndex, 
						INT_4 *aaPresentMass, 
						INT_4 topSeqNum, INT_4 *lastNode, INT_4 lastNodeNum, 
						INT_4 *seqNum, INT_4 maxLastNode, INT_4 minLastNode,
						INT_4 lowSuperNode, INT_4 highSuperNode)
{
	struct Sequence *newSubsequencePtr;
	struct Sequence *currPtr;
	char doIt;		
	BOOLEAN oneEdgeNNode, test;
	INT_4 i, j, k, testValue, subseqNum, z=0, subseqCount;
	INT_4 extNum;
	INT_4 massDiff;
	struct extension	clearExtension;
	struct extension 	*extensionList;
	
	
	clearExtension.gapSize        = 0;
	clearExtension.mass           = 0;
	clearExtension.singleAAFLAG   = 0;
	clearExtension.score          = 0;
	clearExtension.nodeCorrection = 0;

	
	extensionList = (extension *) calloc(MAX_GAPLIST, sizeof(struct extension)); /* Note: memory is zeroed */
	if(extensionList == NULL)
	{
		printf("AddExtensions:  Out of memory");
		exit(1);
	}

	subseqCount = 0;	/*This is used to count the number of subsequences in the old list*/
	newSubsequencePtr = NULL;
	currPtr = subsequencePtr;
	subseqNum = 0;	/*This is used to count the number of subsequences in the new list
					  and is used in the StoreSubsequences function to ensure that 
					  only gParam.topSeqNum subsequences are stored.*/

	while(currPtr != NULL)
	{
		/* Clear the extensionList */		
		for (i = 0; i < MAX_GAPLIST; i++) {
			extensionList[i] = clearExtension;
		}
		

	
		subseqCount++;	
		extNum = 0;
		oneEdgeNNode = TRUE;	/*If it remains TRUE, then no one amino acid extensions found.*/
		for(i = 0; i < gAminoAcidNumber; i++)	/*Find the one amino acid extensions.*/
		{
			if(gGapList[i] != 0)
			{
				testValue = currPtr->nodeValue + gGapList[i];
				doIt = TRUE;
				if(currPtr->nodeValue < lowSuperNode && testValue > highSuperNode)
				{
					doIt = FALSE;
				}
				if(testValue >= gGraphLength)
				{
					doIt = FALSE;	/*you've gone past the graph length*/
				}
				if(sequenceNode[testValue] != 0 && testValue <= maxLastNode && doIt)
				{
					/* Add the single AA extension to the extensionList */
					extensionList[extNum].mass           = gGapList[i];
					extensionList[extNum].gapSize        = 0;
					extensionList[extNum].singleAAFLAG   = 1;
					extensionList[extNum].score          = sequenceNode[testValue];
					extensionList[extNum].nodeCorrection = gNodeCorrection[i];
					extNum++;
					if(extNum >= MAX_GAPLIST)
					{
						printf("AddExtensions:  extNum >= MAX_GAPLIST\n");
						exit(1);
					}
					oneEdgeNNode = FALSE;
				}
			}
		}
		
/*
*	Now find the two amino acid extensions.  Make sure that the two amino acid extension does not
*	include one of the one amino acid extensions.  IE, if I find an edge for Ala, then I cannot
*	use a two amino acid edge of 199, since 128 + 71 = 199.
*/

		for(i = gAminoAcidNumber; i <= gGapListIndex; i++)/*Start at the end of the one aa 
															extensions and move up from there.*/
		{
			testValue = currPtr->nodeValue + gGapList[i];/*This is the test mass (nominal).*/
			doIt = TRUE;	/*check that no superNodes are involved*/
			if(currPtr->nodeValue < lowSuperNode && testValue > highSuperNode)
			{
				doIt = FALSE;
			}
			if(testValue >= gGraphLength)
			{
				doIt = FALSE;	/*you've gone past the graph length*/
			}
			if(doIt && sequenceNode[testValue] != 0 && testValue <= maxLastNode) {
					/* If there is any evidence at that mass,
					   and if the node is less than the peptide mass.*/

				test = TRUE;	/*test is set to true, and if it continues to be true (see below),
							      then this extension is allowed.*/
				
				j = 0;
				while (extensionList[j].singleAAFLAG == 1 && j < extNum) {
					/*Look at the one aa extensions.*/
					massDiff = gGapList[i] - extensionList[j].mass;
					for(k = 0; k < gAminoAcidNumber; k++)	/*k = 0; k < gGapListIndex; k++*/
					{
						if(gGapList[k] != 0)
						{
							if(massDiff <= gGapList[k] + gParam.fragmentErr
								&& massDiff >= gGapList[k] - gParam.fragmentErr)
							{
								test = FALSE;	/*If the mass difference is equal to an amino acid
												mass then its not allowed and test = FALSE.*/
								break;
							}
						}
					}
					
					j++;
				}				
				if(test)	/*If test is still true then go ahead and save these as extensions.*/
				{
					/*Add the 2 AA extension.*/
					extensionList[extNum].mass           = gGapList[i];
					extensionList[extNum].gapSize        = 1;
					extensionList[extNum].singleAAFLAG   = 0;
					extensionList[extNum].nodeCorrection = 0;
					
					/*Examine some special cases and then calculate the score for the 2 AA extension*/
					extensionList[extNum] = Score2aaExtension(extensionList[extNum], 
															  currPtr->nodeValue,
															  oneEdgeNodes,
															  oneEdgeNodesIndex,
															  oneEdgeNNode,
															  sequenceNode[testValue]);					
					
					extNum++;	/*Increment the number of extensions.*/
					if(extNum >= MAX_GAPLIST)
					{
						printf("AddExtensions:  extNum >= MAX_GAPLIST\n");
						exit(1);
					}
				}
			}
		}

/*
*	Now I need to find the best extensions, ie, the top maxExtNum of them and only if these
*	extensions are greater than the product of the highest score and extThresh.
*	 bestExtensions[MAX_GAPLIST], bestExtScore[MAX_GAPLIST], bestExtNum;
*/
		
		if(extNum > 0)
		{
			if(currPtr->peptideLength >= MAX_PEPTIDE_LENGTH)
			{
				printf("LutefiskSubseqMaker:AddExtensions  peptide length exceeds array length");
				exit(1);
			}
			extensionList = SortExtensions( extensionList );

			
					
/*	Store this information in the linked list of Sequence structs.  The values placed in the
*	arrays bestExtensions[], bestExtScore[], and bestExtGapNum[] are used to determine
*	values for the variables peptide[], score, peptideLength, gapNum, and nodeValue.  These
*	variables are passed to a few functions that are used to set up the linked list of structs
*	of type Sequence, which contain the next set of subsequences (newSubsequencePtr*/

			newSubsequencePtr = StoreSubsequences(newSubsequencePtr, extensionList, currPtr,
												  lastNode, lastNodeNum, maxLastNode, minLastNode,
												  aaPresentMass, seqNum, &subseqNum,
												  sequenceNode);

			
		}
		currPtr = currPtr->next;	/*Go to the subsequence.*/
	}
	
/*
*	Hack back the number of subsequences to be less than gParam.topSeqNum.
*/

	/*while(subseqNum > gParam.topSeqNum)
	{
		ClearLowestScore(newSubsequencePtr, &subseqNum, maxLastNode, minLastNode);
	}*/
	
	if(subseqCount > gSubseqNum)
	{
		gSubseqNum = subseqCount;	/*keep track of the most numbers of subsequences used*/
	}
	
	free(extensionList);

	FreeSequenceStructs(subsequencePtr);

	return(newSubsequencePtr);
}

/****************LinkSubsequenceList**********************************************************
*
* 	This function adds a subsequence onto the existing linked list of structs of type Sequence.
*	It adds structs in order of their score fields, so that the first in the list has the
*	highest score and the last in the list has the lowest score.
*
*/

struct Sequence *LinkSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr)
{
	struct Sequence *currPtr, *previousPtr;
	char test = TRUE;
		
	if(firstPtr == NULL)	/*If this is the first struct of the list then do this.*/
		firstPtr = newPtr;
	else
	{
		if(newPtr->score > firstPtr->score)	/*If the struct to be added has the best score then*/
		{
			newPtr->next = firstPtr;
			firstPtr = newPtr;
		}
		else	/*Otherwise, go find the position for the new struct (based on its score field.*/
		{
			previousPtr = firstPtr;
			currPtr = firstPtr->next;
			while(currPtr != NULL)
			{
				if(newPtr->score > currPtr->score)	/*I found the place.*/
				{
					previousPtr->next = newPtr;
					newPtr->next = currPtr;
					test = FALSE;
					break;
				}
				previousPtr = currPtr;
				currPtr = currPtr->next;
			}
			if(test)	/*Ok, I didn't find the place, and I'm at the end of the linked list,
						so this must have the lowest score of all.  I'll stick it on the end.*/
			{
				previousPtr->next = newPtr;
				newPtr->next = NULL;
			}
		}
	}
	
	return(firstPtr);
}

/******************LoadFinalSequenceStruct********************************************
*
* 	LoadStruct puts the nominal extension mass in the peptide[] field, and increments the 
*	value of peptideLength by one.  The fields score and nodeValue are also modified.  Of
*	course, to do this the function finds some memory, and this value is returned as a pointer
*	to a struct of type Sequence (which now contains all of this data).
*	struct Sequence
*	{
*		INT_4 peptide[MAX_PEPTIDE_LENGTH];
*		INT_4 peptideLength;
*		INT_4 score;
*		INT_4 nodeValue;
*		struct Sequence *next;
*	};
*/
struct Sequence *LoadFinalSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
							INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection)
{
	struct Sequence *currPtr;
	INT_4 i;	
	REAL_4 scoreAdjuster;

	currPtr = (struct Sequence *) malloc(sizeof(struct Sequence));
	if(currPtr == NULL)
	{
		printf("LoadFinalSequenceStruct:  Out of mammories");
		exit(1);
	}
	
	scoreAdjuster = nodeValue;				/*Here's the expected average peptide length.*/
	scoreAdjuster = scoreAdjuster / gAvResidueMass;
	if(peptideLength == 0)
	{
		printf("LoadFinalSequenceStruct:  peptideLength == 0\n");
		exit(1);
	}
	scoreAdjuster = scoreAdjuster / peptideLength;	/*Here's the ratio between the average
													expected and actual length.*/
	
	for(i = 0; i < peptideLength; i++)
	{
		currPtr->peptide[i] = peptide[i];
	}	
	currPtr->peptideLength = peptideLength;	
	currPtr->score = score * scoreAdjuster;
	currPtr->gapNum = gapNum;
	currPtr->nodeValue = nodeValue;
	currPtr->nodeCorrection = nodeCorrection;
	currPtr->next = NULL;

	return(currPtr);
}

/******************LoadSequenceStruct********************************************
*
* 	LoadStruct puts the nominal extension mass in the peptide[] field, and increments the 
*	value of peptideLength by one.  The fields score and nodeValue are also modified.  Of
*	course, to do this the function finds some memory, and this value is returned as a pointer
*	to a struct of type Sequence (which now contains all of this data).
*	struct Sequence
*	{
*		INT_4 peptide[MAX_PEPTIDE_LENGTH];
*		INT_4 peptideLength;
*		INT_4 score;
*		INT_4 nodeValue;
*		struct Sequence *next;
*	};
*/
struct Sequence *LoadSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
							INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection)
{
	struct Sequence *currPtr;
	INT_4 i;	

	currPtr = (struct Sequence *)calloc(1, sizeof(struct Sequence));
	if(currPtr == NULL)
	{
		printf("LoadSequenceStruct:  Out of mammories");
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

/**********************************NterminalSubsequences****************************************
*
*	This function sets up the first batch of subsequences.  It starts with a single subsequence
*	containing the N-terminal group (usually hydrogen of mass 1) and then tries to connect with
*	nodes that are either one or two amino acids higher in mass.  It uses the array gapList
*	to determine what is an acceptable mass difference.  A linked list of structs of type
*	Sequence is generated where the first struct in the list has the highest score, and the last
*	struct in the list has the lowest score.  There is an upper limit on the number of 
*	N-terminal extensions - maxExtNum - and the scores for these must be greater than the
*	product of the highest scoring subsequence and the value of extThresh.  If the peptide
*	mass is small enough, then the sequences are stored in a linked list of final sequences,
*	of which there will be an upper limit of finalSeqNum.  This function returns a pointer to 
*	a struct of type Sequence, which is the first element in a linked list of subsequences.
*/

struct Sequence *NterminalSubsequences(SCHAR *sequenceNode, INT_4 maxLastNode, INT_4 lowSuperNode,
										INT_4 highSuperNode)
{
	struct Sequence *subsequencePtr;
	INT_4 i, j, k, m, testValue, nTerminus;
	INT_4 *extensions, *extScore, extNum;
	INT_4 *bestExtensions, *bestExtScore, bestExtNum;
	INT_4 highestExtensionScore;
	INT_4 threshold;
	INT_4 score, *peptide, peptideLength, nodeValue, gapNum;
	INT_4 *threeAA, threeAANum, sum;
	INT_4 sameNum, averageExtension, *sameExtension;
	INT_2 nodeCorrection;
	char nTerminusPossible, duplicateFlag;
	char doIt;
	char sameTest;
	
	sameExtension = (int *) malloc(gParam.fragmentErr * 20 * sizeof(INT_4));
	if(sameExtension == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}
	
	threeAA = (int *) malloc(gAminoAcidNumber*gAminoAcidNumber*gAminoAcidNumber*sizeof(INT_4));
	if(threeAA == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}
	extensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(extensions == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}
	extScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(extScore == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}
	bestExtensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(bestExtensions == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}
	bestExtScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4 ));
	if(bestExtScore == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}

	peptide = (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4 ));
	if(peptide == NULL)
	{
		printf("NterminalSubsequences:  Out of memory.");
		exit(1);
	}

/*	Fill in the masses for three amino acids.*/

	threeAANum = 0;
	for(i = 0; i < gAminoAcidNumber; i++)	/*Fill in the masses of the 3 AA extensions.*/
	{
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			for(k = 0; k < gAminoAcidNumber; k++)
			{
				if(gGapList[i] != 0 && gGapList[j] != 0 && gGapList[k] != 0)
				{
					sum = gGapList[i] + gGapList[j] + gGapList[k];
					duplicateFlag = FALSE;
					for(m = 0; m < threeAANum; m++)
					{
						if(threeAA[m] == sum)
						{
							/* We already have this mass in threeAA so don't add it to the list. */
							duplicateFlag = TRUE;
							break;
						}
					}
					if(duplicateFlag == FALSE)
					{
						for(m = 0; m <= gGapListIndex; m++)
						{
							if(gGapList[m] == sum)
							{
							/* We already have this mass  so don't add it to the list. */
							duplicateFlag = TRUE;
							break;
							}
						}
					}
					if(!duplicateFlag)
					{
						threeAA[threeAANum] = sum;
						threeAANum++;
					}
				}
			}
		}
	}

	
	extNum = 0;
	subsequencePtr = NULL;
	nTerminus = gParam.modifiedNTerm + 0.5;
	


/*	Find the one, two, and three amino acid jumps from gGapList and threeAA.*/

	for(i = nTerminus + gMonoMass_x100[G]; i < gMonoMass_x100[W] * 3; i++)	/*step thru each node*/
	{
		if(sequenceNode[i] != 0)	/*ignore the nodes w/ zero evidence*/
		{
			nTerminusPossible = FALSE;	/*start assuming that this is not an extension*/
/*	Check for one and two amino acid extensions using the gGapList array.*/
			for(j = 0; j <= gGapListIndex; j++)	
			{
				if(gGapList[j] != 0)
				{
					if(i - nTerminus == gGapList[j])
					{
						nTerminusPossible = TRUE;	/*its a possible extension*/
						break;
					}
				}
			}
/*	Now check for three amino acid extensions if it wasn't a one or two aa extension.*/
			if(nTerminusPossible == FALSE)
			{
				for(j = 0; j < threeAANum; j++)
				{
					if(threeAA[j] != 0)
					{
					 	if(i - nTerminus == threeAA[j])
					 	{
					 		nTerminusPossible = TRUE;
					 		break;
					 	}
					}
				}
			}
			doIt = TRUE;	/*check for superNodes when sequencetag specified*/
			if(nTerminus < lowSuperNode && i > highSuperNode)
			{
				doIt = FALSE;
			}
			if(nTerminusPossible && i <= maxLastNode && doIt)	/*save as an extension?*/
			{
				extensions[extNum] = i - nTerminus;
				extScore[extNum] = sequenceNode[i];
				extNum++;
			}
		}
	}
	
/*	Find extensions that are 1 node unit apart, and consolidate them.*/

	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] != 0)
		{
			sameNum = 0;
			averageExtension = extensions[i];
			for(j = 0; j < extNum; j++)
			{
				if(j != i && extScore[j] != 0)
				{
					sameTest = FALSE;
					for(k = 0; k < sameNum; k++)
					{
						if(extensions[sameExtension[k]] - extensions[j] == 1 ||
							extensions[j] - extensions[sameExtension[k]] == 1)
						{
							sameTest = TRUE;
						}
					}
					if(extensions[i] - extensions[j] == 1 ||
						extensions[j] - extensions[i] == 1 || sameTest)
					{
						sameExtension[sameNum] = j;
						averageExtension += extensions[j];
						sameNum++;
					}
				}
			}
			if(sameNum != 0)
			{
				averageExtension = ((float)averageExtension / (sameNum + 1)) + 0.5;	/*count the i extension
																			and round the value*/
				extensions[i] = averageExtension;	
				for(j = 0; j < sameNum; j++)
				{
					extScore[sameExtension[j]] = 0;
				}
			}
		}
	}			
	
/*	Get rid of the extensions that are within gParam.fragmentErr of each other.*/
	
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] != 0)
		{
			for(j = 0; j < extNum; j++)
			{
				if(extScore[j] != 0 && i != j)
				{
					if(extensions[i] <= extensions[j] + gParam.fragmentErr
						&& extensions[i] >= extensions[j] - gParam.fragmentErr)
					{
						if(extScore[i] >= extScore[j])
						{
							extScore[j] = 0;
						}
						else
						{
							extScore[i] = 0;
						}

					}
				}
			}
		}
	}
		
	
/*	If two extensions differ by the mass of an amino acid, then the higher mass one's 
	intensity is assigned a zero so that it is removed in the section below.*/
	
	if(extNum > 0)
	{
		for(i = 0; i < extNum; i++)
		{
			for(j = i + 1; j < extNum; j++)
			{
				testValue = extensions[j] - extensions[i];
				for(k = 0; k < gAminoAcidNumber; k++)
				{
					if(gGapList[k] != 0)
					{
						if(testValue <= gGapList[k] + gParam.fragmentErr &&
							testValue >= gGapList[k] - gParam.fragmentErr)
						{
							extScore[j] = 0;
						}
					}
				}
			}
		}
	}
	
/*
*	Now I need to find the best extensions, ie, the top maxExtNum of them and only if these
*	extensions are greater than the product of the highest score and extThresh.
*	 bestExtensions[MAX_GAPLIST], bestExtScore[MAX_GAPLIST], bestExtNum;
*/

	highestExtensionScore = extScore[0];	/*Find the highest extension score.*/
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] > highestExtensionScore)
		{
			highestExtensionScore = extScore[i];
		}
	}
	
	threshold = highestExtensionScore * gParam.extThresh; /*Set the extension score threshold.*/
/*	
*	If a peptide tag has been entered and if the highestExtensionScore is over 100, then 
*	that means that the highest scoring extension is to a superNode.  In this case, don't
*	use an extension score threshold.
*/
	if(gParam.tagNMass != 0 && gParam.tagCMass != 0)
	{
		if(highestExtensionScore > 100)
		{
			threshold = 0;
		}
	}


	bestExtNum = 0;
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] >= threshold)
		{
			bestExtensions[bestExtNum] = extensions[i];
			bestExtScore[bestExtNum] = extScore[i];
			bestExtNum++;
		}
	}
	
/*	Make sure there are enough subsequences allowed.*/
	if(bestExtNum > gParam.topSeqNum)
	{
		bestExtNum = gParam.topSeqNum;
	}
			
/*	Store this information in the linked list of Sequence structs.*/

	for(i = 0; i < bestExtNum; i++)
	{
		score = bestExtScore[i];
		peptide[0] = bestExtensions[i];
		peptideLength = 1;
		gapNum = 0;
		nodeValue = nTerminus + bestExtensions[i];
		j = gParam.modifiedNTerm * 10 + 0.5;
		k = gParam.modifiedNTerm + 0.5;
		k = k * 10;
		nodeCorrection = j - k;

		subsequencePtr = LinkSubsequenceList(subsequencePtr, LoadSequenceStruct(peptide, 
																	peptideLength, 
																	score, nodeValue, gapNum,
																	nodeCorrection));
	}
	
	free(extensions);
	free(extScore);
	free(bestExtensions);
	free(bestExtScore);
	free(peptide);
	free(sameExtension);
	free(threeAA);
	
	return(subsequencePtr);
}


/**********************************SubsequenceMaker********************************************
*
*	This function uses a subsequencing approach to derive a list of completed peptide sequences.
*	The array sequenceNode contains the information used to build the subsequences.  The indexing
*	of this array corresponds to nominal masses of b ions, and the information contained 
*	relates to the probability that a cleavage was actually present at that nominal mass.  The
*	array oneEdgeNodes lists the nodes that could not be connected to any other node N-terminal
*	to it, and oneEdgeNodesIndex is the number of such nodes listed in the array.  The arrays
*	aaPresent and aaAbsent contain the amino acids (single letter code) that are either present
*	or absent in the peptide.  The value of cysMW allows for alteration of the residue mass of
*	cysteine, which can be alkylated w/ various reagents.  The extThresh is the fractional
*	value of the highest ranked extension to a particular subsequence that can be used in 
*	the formation of the next list of subsequences.  The value of maxExtNum is the maximum
*	number of extensions per subsequence that is allowed.  The finalSeqNum is the upper limit
*	on the number of completed sequences that will be stored and scored later.  The topSeqNum
*	is the upper limit on the number of subsequences that are to be allowed.  If I run out
*	of space, I may make these last two parameters changeable by the program - ie, if there
*	is no space left then it will automatically purge 100 or so sequences and reset these
*	upper limits to 100 less than before.  The maxGapNum is the maximum number of gaps that 
*	are allowed - its either zero or one.
*
*/

struct Sequence *SubsequenceMaker(INT_4 *oneEdgeNodes, INT_4 oneEdgeNodesIndex, 
						SCHAR *sequenceNode)
{
	struct Sequence *subsequencePtr;
	INT_4 *lastNode, lastNodeNum, aaPresentMass[AMINO_ACID_NUMBER];
	INT_4 maxLastNodeNum;	/*max num of lastnodes*/
	INT_4 maxLastNode, minLastNode;		/*the highest and lowest last node value*/
	INT_4 i, j, seqNum, finalSeqNum, topSeqNum, correctPeptideLength;
	INT_4 highSuperNode, lowSuperNode;	/*used when a specific sequence tag is to be used*/
	INT_4 halfAsManySubsequences, quarterAsManySubsequences;
	char test;
	
	gFinalSequencePtr = NULL;
	subsequencePtr = NULL;
	seqNum = 0;
	gSubseqNum = 0;
	finalSeqNum = gParam.finalSeqNum;
	topSeqNum = gParam.topSeqNum;
	halfAsManySubsequences		= gParam.topSeqNum / 2;
	quarterAsManySubsequences	= gParam.topSeqNum / 4;


/*	Determine the values of highSuperNode and lowSuperNode*/
	highSuperNode = gGraphLength -1;
	lowSuperNode = 0;
	
	if(gParam.tagSequence[0] != '*')
	{
		for(i = gGraphLength - 1; i >= 0; i--)
		{
			if(sequenceNode[i] == -1)
			{
				highSuperNode = i;
				break;
			}
		}
		
		for(i = 0; i < gGraphLength; i++)
		{
			if(sequenceNode[i] == -1)
			{
				lowSuperNode = i;
				break;
			}
		}
		for(i = 0; i < gGraphLength; i++)
		{
			if(sequenceNode[i] == -1)
			{
				sequenceNode[i] = 127;
			}
		}
	}
/*	Calculated the maximum number of last nodes (this was a source of a memory over-run.*/

	if(gParam.fragmentErr > gParam.peptideErr)
	{
		maxLastNodeNum = 4 * gParam.fragmentErr;
	}
	else
	{
		maxLastNodeNum = 4 * gParam.peptideErr;
	}
	
	lastNode = (int *) malloc(maxLastNodeNum * sizeof(INT_4 ));	/*Will contain C-terminal evidence.*/
	if(lastNode == NULL)
	{
		printf("SubsequenceMaker:  Out of memory");
		exit(1);
	}

	
/*	Determine values for aaPresentMass (nominal masses for the single letter code found in
	aaPresent.*/
	if(gParam.aaPresent[0] == '*')
	{
		aaPresentMass[0] = -1;
	}
	else
	{
		i = 0;
		while(gParam.aaPresent[i] != 0)
		{
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gParam.aaPresent[i] == gSingAA[j])
				{
					if(gGapList[j] != 0)
					{
						aaPresentMass[i] = gGapList[j];
					}
					break;
				}
			}
			i++;
		}
		aaPresentMass[i] = -1;
	}
/*
*	Find the C-terminal nodes, so that I'll know when I'm finished.
*/


	lastNodeNum = 0;
	i = gGraphLength - 1;
	test = TRUE;
	while(i > 0 && test)
	{
		if(sequenceNode[i] > 0)
		{
			test = FALSE;
			maxLastNode = i;
			while(sequenceNode[i] > 0)
			{
				lastNode[lastNodeNum] = i;
				minLastNode = i;
				lastNodeNum++;
				if(lastNodeNum > maxLastNodeNum)
				{
					printf("lastNodeNum exceeds maxLastNodeNum.");
					exit(1);
				}
				i--;
			}
		}
		i--;
	}
	
	if(gParam.proteolysis == 'T')
	{
		gAA1Max = maxLastNode - gMonoMass_x100[R];
		gAA1Min = minLastNode - gMonoMass_x100[R];
		gAA1 = gMonoMass_x100[R];
		gAA2Max = maxLastNode - gMonoMass_x100[K];
		gAA2Min = minLastNode - gMonoMass_x100[K];
		gAA2 = gMonoMass_x100[K];
	}
	else if(gParam.proteolysis == 'K')
	{
		gAA1Max = maxLastNode - gMonoMass_x100[K];
		gAA1Min = minLastNode - gMonoMass_x100[K];
		gAA1 = gMonoMass_x100[K];
		gAA2Max = 0;
		gAA2Min = 0;
		gAA2 = 0;
	}
	else if(gParam.proteolysis == 'E')
	{
		gAA1Max = maxLastNode - gMonoMass_x100[E];
		gAA1Min = minLastNode - gMonoMass_x100[E];
		gAA1 = gMonoMass_x100[E];
		gAA2Max = maxLastNode - gMonoMass_x100[D];
		gAA2Min = minLastNode - gMonoMass_x100[D];
		gAA2 = gMonoMass_x100[D];
	}

	
			
/*
*	The function NterminalSubsequences starts at the N-terminal node (usually a value of one),
*	and jumps by one amino acid, or by two amino acids to generate the first set of subsequences.
*	This linked list of subsequences is passed as a pointer to the first element in the array.
*/
	
	if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q')
	{
		subsequencePtr = NterminalSubsequences(sequenceNode, maxLastNode, lowSuperNode,
												highSuperNode);
	}
	if(gParam.fragmentPattern == 'L')
	{
		subsequencePtr = LCQNterminalSubsequences(sequenceNode, maxLastNode, lowSuperNode,
												highSuperNode);
	}
	
	if(gCheckItOut)	/*If I want to see how well my subsequencing is going gCheckItOut is TRUE.*/
	{
		correctPeptideLength = 1;
		amIHere(correctPeptideLength, subsequencePtr);
	}
							
/*
*	The function AddExtensions primarily adds one amino acid at a time to the existing list
*	of subsequences.  If a subsequence cannot be extended, then a two amino acid gap is jumped
*	to one of the one-edged nodes from the C-terminus.  No penalty is extracted as a consequence
*	of jumping, since its only allowed between two nodes that cannot be extended.  Other
*	two amino acid jumps are penalized so that the extension score is reduced by factors 
*	NODE_NODE_PENALTY, NODE_EDGE_PENALTY, and EDGE_EDGE_PENALTY as defined in lutefisk.h.
*	Once there are no more nodes remaining, then the function returns a NULL value.
*/
	
	while(subsequencePtr != NULL)
	{
		if((clock() - gParam.startTicks)/ CLOCKS_PER_SEC > 30)
		{
			gParam.topSeqNum = halfAsManySubsequences;	/*its taking too long, so speed it up*/
		}
		if((clock() - gParam.startTicks)/ CLOCKS_PER_SEC > 60)
		{
			gParam.topSeqNum = quarterAsManySubsequences;	/*this is really taking too long*/
		}
		
		subsequencePtr = AddExtensions(subsequencePtr, sequenceNode, 
									oneEdgeNodes, oneEdgeNodesIndex, aaPresentMass, 
									topSeqNum, lastNode, lastNodeNum, &seqNum,
									maxLastNode, minLastNode, lowSuperNode, highSuperNode);
	
		if(gCheckItOut)
		{
			correctPeptideLength++;
			amIHere(correctPeptideLength, subsequencePtr);
		}
	}
	
	if(gCheckItOut)
	{
		correctPeptideLength = 14;
		amIHere(correctPeptideLength, gFinalSequencePtr);
	}
	if(gParam.fMonitor && gCorrectMass)
	{
		printf("Subsequencing is finished.\n");
		printf("Max subsequences: %d \n", gSubseqNum);
	}

	free(lastNode);
	
	return(gFinalSequencePtr);
}
