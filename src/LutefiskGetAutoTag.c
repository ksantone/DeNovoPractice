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

#include <stdio.h>
#include <stdlib.h>
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"
#include "ListRoutines.h"


/*Globals for this file only.*/
tsequenceList *TagSeqList;
tsequenceList *TagSubseqList;

tMSDataList   *TagMassList;

/************************************* SortExtension **************************************
*
*
*
*/

struct extension *SortExtension(struct extension *inExtensionList)
{
	INT_4 i;
	INT_4 extensionCount, outIndex;
	INT_4 highestExtensionScore, threshold;
	REAL_4 extensionThreshold = 0.05;	/*extensions must have score greater than this*/
	
	
	struct extension *outExtensionList;
	
	outExtensionList = (extension *) calloc(MAX_GAPLIST, sizeof(struct extension));
	if(outExtensionList == NULL)
	{
		printf("SortExtension:  Out of memory.");
		exit(1);
	}
	
	/*Find the highest extension score.*/
	highestExtensionScore = inExtensionList[0].score; 
	extensionCount = 0;	
	i = 0;
	while(inExtensionList[i].mass > 0 && i < MAX_GAPLIST)
	{
		if(inExtensionList[i].score > highestExtensionScore)
		{
			highestExtensionScore = inExtensionList[i].score;
		}
		i++;
		extensionCount++;
	}
	
	/*Set the extension score threshold.*/
	threshold = highestExtensionScore * extensionThreshold; 	
	
	/*If the number of extensions found are less than the specified
	  maximum number of extensions per subsequence, then do this.*/
	if(extensionCount <= gParam.maxExtNum)	
	{
		outIndex = 0;
		for(i = 0; i < extensionCount; i++)
		{
			if(inExtensionList[i].score >= threshold && i < gParam.topSeqNum)
			{
				outExtensionList[outIndex] = inExtensionList[i];
				outIndex++;
				if(outIndex >= MAX_GAPLIST)
				{
					printf("SortExtension:  outIndex >= MAX_GAPLIST\n");
					exit(1);
				}
			}
		}
	}
	if(extensionCount > gParam.maxExtNum)
	{
		/*There are too many extensions; throw out the worst ones.*/
		
		/*First sort the extensions by whether they are singleAA or not and then by score.*/
		qsort(inExtensionList, extensionCount, sizeof(struct extension), ExtensionsSortDescend);
		
		outIndex = 0;
		i = 0;
		while(i < extensionCount && outIndex < gParam.maxExtNum && outIndex < gParam.topSeqNum)
		{
			if(inExtensionList[i].score >= threshold)
			{
				outExtensionList[outIndex] = inExtensionList[i];
				outIndex++;
				if(outIndex >= MAX_GAPLIST)
				{
					printf("SortExtensions:  outIndex >= MAX_GAPLIST\n");
					exit(1);
				}
			}
			i++;
		}
	}
	
	free(inExtensionList);
	
	return(outExtensionList);
}	


/***************************RemoveNeutralLosses****************************************
*
*	Anything that is 17 or 18 Da less than another ion, where the neutral loss ion is
*	less than half the intensity of the higher mass ion is discarded.
*/
void RemoveNeutralLosses(INT_4 charge)
{
	REAL_4 massDiff;
	char test;
	tMSData *tagMassPtr, *nextMassPtr;
		
	if(TagMassList->numObjects < 2)
	{
		return;
	}
	
	tagMassPtr = &TagMassList->mass[0];
	while(tagMassPtr < TagMassList->mass + (TagMassList->numObjects) - 1)
	{
		nextMassPtr = tagMassPtr + 1;
		while(nextMassPtr < TagMassList->mass + TagMassList->numObjects)
		{
			test = FALSE;
			massDiff = nextMassPtr->mOverZ - tagMassPtr->mOverZ;
			massDiff = massDiff * charge;	/*deal w/ charge > 1*/
			if(massDiff <= gWater + gParam.fragmentErr && 
				massDiff >= gAmmonia - gParam.fragmentErr)
			{
				if((2 * tagMassPtr->intensity) < nextMassPtr->intensity)
				{
					RemoveFromList(tagMassPtr - TagMassList->mass, TagMassList);
					break;
				}
			}
			nextMassPtr++;
		}
		tagMassPtr++;
	}
	return;
}

/***************************FilterTagMasses********************************************
*
*	Get rid of ions that are closer than Gly (remove the lower intensity of the pair).
*
*/
void FilterTagMasses(INT_4 charge)
{
	REAL_4 massDiff, minDiff;
	
	tMSData *tagMassPtr;
	
	
	minDiff = (gGapList[G] / charge) - 2 * gParam.fragmentErr;	/*2x err for a little slop*/
	
	if(TagMassList->numObjects < 2) 
	{
		return;
	}

/*If things are too close, then zero out the intensity of the lowest ones*/	
	tagMassPtr = &TagMassList->mass[1];	/*start w/ the second tagMass*/
	while(tagMassPtr < TagMassList->mass + TagMassList->numObjects) 
	{
		massDiff = tagMassPtr->mOverZ - (tagMassPtr - 1)->mOverZ;
		if(massDiff < minDiff)
		{
			if((2 * (tagMassPtr - 1)->intensity) < tagMassPtr->intensity)
			{	/*just remove things that are more than 2x diff in intensity*/
				RemoveFromList(tagMassPtr - 1 - TagMassList->mass, TagMassList); 
				tagMassPtr--;
			}
			else if((2 * tagMassPtr->intensity) < (tagMassPtr - 1)->intensity)
			{
				RemoveFromList(tagMassPtr - TagMassList->mass, TagMassList); 
				tagMassPtr--;
			}
		}
		
		tagMassPtr++;
	}
}
/**********************************NterminalTags****************************************
*
*	This function sets up the first batch of subsequences.  It starts with a single subsequence
*	containing the N-terminal group (usually hydrogen of mass 1) and then tries to connect with
*	nodes that are either one or two amino acids higher in mass.  It uses the array gapList
*	to determine what is an acceptable mass difference.  A linked list of structs of type
*	Sequence is generated where the first struct in the list has the highest score, and the last
*	struct in the list has the lowest score.  This function returns a pointer to 
*	a struct of type Sequence, which is the first element in a linked list of subsequences.
*/

INT_4 NterminalTags(char *tagNode, INT_4 *tagNodeIntensity)
{
	struct Sequence *subsequencePtr;
	tsequence subsequenceToAdd;
	INT_4 i, j, k, testValue, nTerminus;
	INT_4 *extensions, *extScore, extNum;
	INT_4 *bestExtensions, *bestExtScore, bestExtNum;
	INT_4 highestExtensionScore;
	INT_4 *singleAAExtension, singleAAExtNum, threshold;
	INT_4 *peptide;
	
	extensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(extensions == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}
	
	extScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(extScore == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}

	bestExtensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(bestExtensions == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}

	bestExtScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(bestExtScore == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}

	singleAAExtension = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(singleAAExtension == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}

	peptide = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(peptide == NULL)
	{
		printf("NterminalTags:  Out of memory");
		exit(1);
	}

	
	extNum = 0;
	singleAAExtNum = 0;
	subsequencePtr = NULL;
	
	/*Figure out what the N-terminal node is.*/
	nTerminus = gParam.modifiedNTerm + 0.5;
	
	for(i = 0; i < gAminoAcidNumber; i++)	/*Find the one amino acid extensions.*/
	{
		if(gGapList[i] != 0)
		{
			testValue = nTerminus + gGapList[i];
			if(tagNode[testValue] != 0)
			{
				extensions[extNum] = gGapList[i];
				singleAAExtension[extNum] = gGapList[i];
				extScore[extNum] = tagNodeIntensity[testValue];
				extNum++;
				singleAAExtNum++;
			}
		}
	}
	
/*
	Now find the two amino acid extensions.
*/

	for(i = gAminoAcidNumber; i <= gGapListIndex; i++)/*Start at the end of the one aa 
															extensions and move up from there.*/
	{
		testValue = nTerminus + gGapList[i];/*This is the test mass (nominal).*/
		if(tagNode[testValue] != 0)	/*If there is any evidence at that mass.*/
		{
			extensions[extNum] = gGapList[i];
			extScore[extNum] = tagNodeIntensity[testValue];
			extNum++;
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
						if(testValue == gGapList[k])
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
	
	threshold = highestExtensionScore * 0.01; /*Set the extension score threshold.*/

	bestExtNum = 0;
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] > threshold)
		{
			bestExtensions[bestExtNum] = extensions[i];
			bestExtScore[bestExtNum] = extScore[i];
			bestExtNum++;
		}
	}
			
/*	Store this information in the linked list of Sequence structs.*/

	for(i = 0; i < bestExtNum; i++)
	{
		subsequenceToAdd.score = bestExtScore[i];
		subsequenceToAdd.peptide[0] = bestExtensions[i];
		subsequenceToAdd.peptideLength = 1;
		subsequenceToAdd.gapNum = 0;
		subsequenceToAdd.nodeValue = nTerminus + bestExtensions[i];
		if(!AddToList(&subsequenceToAdd, TagSubseqList)) 
		{
			free(extensions);
			free(extScore);
			free(bestExtensions);
			free(bestExtScore);
			free(singleAAExtension);
			free(peptide);
			return(extNum);
		}
	}
	
	free(extensions);
	free(extScore);
	free(bestExtensions);
	free(bestExtScore);
	free(singleAAExtension);
	free(peptide);
	
	return(extNum);
}
/**********************************AlternateNterminalTags****************************************
*
*	This function sets up the first batch of subsequences.  It starts with a single subsequence
*	containing the N-terminal group (usually hydrogen of mass 1) and then tries to connect with
*	nodes that are either one or two amino acids higher in mass.  It uses the array gapList
*	to determine what is an acceptable mass difference.  A linked list of structs of type
*	Sequence is generated where the first struct in the list has the highest score, and the last
*	struct in the list has the lowest score.  This function returns the number of extensions.
*/

INT_4 AlternateNterminalTags(char *tagNode, INT_4 *tagNodeIntensity)
{
	tsequence subsequenceToAdd;
	INT_4 i, j, k, m, testValue, nTerminus;
	INT_4 *extensions, *extScore, extNum;
	INT_4 *bestExtensions, *bestExtScore, bestExtNum;
	INT_4 highestExtensionScore;
	INT_4 threshold;
	/*INT_4 *peptide;*/
	INT_4 *threeAA, threeAANum, sum;
	INT_4 sameNum, averageExtension, *sameExtension;
	char nTerminusPossible, duplicateFlag;
	char sameTest;
	
	sameExtension = (int *) malloc(gParam.fragmentErr * 10 * sizeof(INT_4));
	if(sameExtension == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}
	
	threeAA = (int *) malloc(gAminoAcidNumber*gAminoAcidNumber*gAminoAcidNumber*sizeof(INT_4));
	if(threeAA == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}

	extensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(extensions == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}

	extScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(extScore == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}

	bestExtensions = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(bestExtensions == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}

	bestExtScore = (int *) malloc(MAX_GAPLIST * sizeof(INT_4));
	if(bestExtScore == NULL)
	{
		printf("AlternateNterminalTags:  Out of memory.");
		exit(1);
	}

	/*peptide = malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4 ));
	if(peptide == NULL)
	{
		printf("Out of memory");
		exit(1);
	}*/

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
							duplicateFlag = true;
							break;
							}
						}
					}
					if(!duplicateFlag)
					{
						threeAA[threeAANum] = sum;
						threeAANum++;
						if(threeAANum >= gAminoAcidNumber*gAminoAcidNumber*gAminoAcidNumber)
						{
							printf("AlternateNTerminalTags:  threeAANum >= gAminoAcidNumber*gAminoAcidNumber*gAminoAcidNumber\n");
							exit(1);
						}
					}
				}
			}
		}
	}
	
/*	Find the n-terminus.*/
	extNum = 0;
	/*Figure out what the N-terminal node is.*/
	nTerminus = gParam.modifiedNTerm + 0.5;

/*	Find the one, two, and three amino acid jumps from gGapList and threeAA.*/

	for(i = nTerminus + gMonoMass_x100[G]; i < gMonoMass_x100[W] * 3; i++)	/*step thru each node*/
	{
		if(tagNode[i] != 0 && i < gGraphLength)	/*ignore the nodes w/ zero evidence*/
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
					
			if(nTerminusPossible)	/*save as an extension?*/
			{
				extensions[extNum] = i - nTerminus;
				extScore[extNum] = tagNodeIntensity[i];
				extNum++;
				if(extNum >= MAX_GAPLIST)
				{
					printf("AlternateNTerminalTags:  extNum >= MAX_GAPLIST\n");
					exit(1);
				}
			}
		}
	}
	
/*	If there are no extensions that match combinations of three amino acids, then look for anything*/
	for(i = nTerminus + gMonoMass_x100[G]; i < gMonoMass_x100[W] * 3; i++)	/*step thru each node*/
	{
		if(tagNode[i] != 0 && i < gGraphLength)	/*ignore the nodes w/ zero evidence*/
		{
			extensions[extNum] = i - nTerminus;
			extScore[extNum] = tagNodeIntensity[i];
			extNum++;
			if(extNum == MAX_GAPLIST - 1)
				break;
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
						if(sameNum >= gParam.fragmentErr * 10)
						{
							printf("AlternateNTerminalTags:  sameNum >= gParam.fragmentErr * 10\n");
							exit(1);
						}
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
			if(extScore[i] != 0)
			{
				for(j = i + 1; j < extNum; j++)
				{
					if(extScore[j] != 0)
					{
						testValue = extensions[j] - extensions[i];
						for(k = 0; k < gAminoAcidNumber; k++)
						{
							if(gGapList[k] != 0)
							{
								if(testValue <= gGapList[k] + gParam.fragmentErr &&
									testValue >= gGapList[k] - gParam.fragmentErr)
								{
									extScore[j] = -1;
								}
							}
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

	if(extNum >= MAX_GAPLIST)
	{
		printf("AlternateNTerminalTags:  extNum >= MAX_GAPLIST\n");
		exit(1);
	}
	highestExtensionScore = extScore[0];	/*Find the highest extension score.*/
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] > highestExtensionScore)
		{
			highestExtensionScore = extScore[i];
		}
	}
	
	threshold = highestExtensionScore * 0.01; /*Set the extension score threshold.*/

	bestExtNum = 0;
	for(i = 0; i < extNum; i++)
	{
		if(extScore[i] > threshold)
		{
			bestExtensions[bestExtNum] = extensions[i];
			bestExtScore[bestExtNum] = extScore[i];
			bestExtNum++;
			if(bestExtNum >= MAX_GAPLIST)
			{
				printf("AlternateNTerminalTags:  bestExtNum >= MAX_GAPLIST\n");
				exit(1);
			}
		}
	}
			
/*	Store this information in the linked list of Sequence structs.*/

	for(i = 0; i < bestExtNum; i++)
	{
		subsequenceToAdd.score = bestExtScore[i];
		subsequenceToAdd.peptide[0] = bestExtensions[i];
		subsequenceToAdd.peptideLength = 1;
		subsequenceToAdd.gapNum = 0;
		subsequenceToAdd.nodeValue = nTerminus + bestExtensions[i];
		if(!AddToList(&subsequenceToAdd, TagSubseqList)) 
		{
			free(extensions);
			free(extScore);
			free(bestExtensions);
			free(bestExtScore);
			free(sameExtension);
			free(threeAA);
			/*free(peptide);*/
			return(extNum);
		}
	}
	
/*	Free the arrays*/
	free(extensions);
	free(extScore);
	free(bestExtensions);
	free(bestExtScore);
	/*free(peptide);*/
	free(sameExtension);
	free(threeAA);
	
	
	return(extNum);
}
/**********************************TagMaker********************************************
*
*	This function uses a subsequencing approach to derive a list of incomplete peptide sequences.
*	The array tagNode contains the information used to build the subsequences.  The indexing
*	of this array corresponds to nominal masses of b ions, and the information contained 
*	relates to the probability that a cleavage was actually present at that nominal mass.  
*	The array tagNodeIntensity contains the ion intensity of the original ions.  
*
*/

void TagMaker(char *tagNode, INT_4 *tagNodeIntensity, INT_4 totalIonIntensity)
{

	INT_4 topSeqNum;
	INT_4 extensions;
	

	topSeqNum = gParam.topSeqNum;
				
/*
*	The function NterminalSubsequences starts at the N-terminal node (usually a value of one),
*	and jumps by one amino acid, or by two amino acids to generate the first set of subsequences.
*	This linked list of subsequences is passed as a pointer to the first element in the array.
*/

	extensions = NterminalTags(tagNode, tagNodeIntensity);
							
/*
*	The function TagExtensions adds one amino acid at a time to the existing list
*	of subsequences.  If a subsequence cannot be extended, then a score is assigned
*	based on intensity.  If the intensity score exceeds a threshold, then the subsequence
*	is stored as a completed tag. 
*	Once there are no more nodes remaining, then the function returns a NULL value.
*/
							
	while(extensions)
	{
		extensions = TagExtensions(tagNode, topSeqNum, tagNodeIntensity, totalIonIntensity);
	}
	
}
/**********************************AlternateTagMaker********************************************
*
*	This function uses a subsequencing approach to derive a list of incomplete peptide sequences.
*	The array tagNode contains the information used to build the subsequences.  The indexing
*	of this array corresponds to nominal masses of b ions, and the information contained 
*	relates to the probability that a cleavage was actually present at that nominal mass.  
*	The array tagNodeIntensity contains the ion intensity of the original ions.  This differs
*	from TagMaker in that the N-terminal amino acids are not singles or REAL_8s; this allows
*	for a larger jump at the N-terminus of the peptide.
*
*/

void AlternateTagMaker(char *tagNode, INT_4 *tagNodeIntensity, INT_4 totalIonIntensity)
{
	INT_4 topSeqNum;
	INT_4 extensions;


	topSeqNum = gParam.topSeqNum;
					
/*
*	The function NterminalSubsequences starts at the N-terminal node (usually a value of one),
*	and jumps to any value between two and three tryptophan masses.
*	This linked list of subsequences is passed as a pointer to the first element in the array.
*/

	extensions = AlternateNterminalTags(tagNode, tagNodeIntensity);
							
/*
*	The function TagExtensions adds one amino acid at a time to the existing list
*	of subsequences.  If a subsequence cannot be extended, then a score is assigned
*	based on intensity.  If the intensity score exceeds a threshold, then the subsequence
*	is stored as a completed tag. 
*	Once there are no more nodes remaining, then the function returns a NULL value.
*/
							
	while(extensions)
	{
		extensions = TagExtensions(tagNode, topSeqNum, tagNodeIntensity, totalIonIntensity);
	}
	
}



/***************************FindTagB2Ions***********************************************
*
*
*
*/

void FindTagB2Ions(struct MSData *firstMassPtr, INT_4 *totalIntensity, 
					char *tagNode, INT_4 *tagNodeIntensity)
{
	INT_4 ionNum;
	struct MSData *currPtr;
	REAL_4 ionMass[258], massDiff, error, bMass;
	INT_4 ionInt[258], i, j, bMassMax, bMassMin, k;
	
/*	
	Since I'm comparing the mass difference of two adjacent ion m/z's to the mass of CO,
	I am guessing that the error will always be fairly tight; at least better than the
	standard m/z fragment error of 0.5 Da.
*/
	if(gParam.fragmentErr > 0.25)
	{
		error = 0.25;
	}
	else
	{
		error = gParam.fragmentErr;
	}
	
/*	Putting the ions into an array makes it easier to look for ions differring by 28 Da.*/
	ionNum = 0;
	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > (115.0 - gParam.fragmentErr) &&
			currPtr->mOverZ < (372.2 + gParam.fragmentErr))
		{
			ionMass[ionNum] = currPtr->mOverZ;
			ionInt[ionNum] = currPtr->intensity;
			ionNum++;
			if(ionNum >= 258)
			{
				printf("The ionNum is greater than 258.");
				exit(1);
			}
		}
		currPtr = currPtr->next;
	}
	
	for(i = 1; i < ionNum; i++)
	{
		j = i-1;
		while(j > 0)
		{
			massDiff = ionMass[i] - ionMass[j];
			if(massDiff >= 28 - error && massDiff <= 28 + error)
			{
				if(ionInt[j] < ionInt[i] * 0.1 ||
					ionInt[j] > ionInt[i] * 3)
					{
						break;	/*if the a/b pair intensity doesn't look right, 
								then toss it out.*/
					}
				totalIntensity += ionInt[i];	/*Boost the total ion intensity.*/
				bMass = ionMass[i];
/*Find highest mass node.*/
				bMass = bMass + gParam.fragmentErr;
				bMassMax = bMass;		/*Truncate w/o rounding up.*/
				tagNode[bMassMax] = 1;
				tagNodeIntensity[bMassMax] += ionInt[i];
/*Next I'll find the lowest mass node.*/
				bMass = bMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				bMassMin = bMass;	/*Truncate after rounding up.*/
				if(bMassMin != bMassMax)
				{
					tagNode[bMassMin] = 1;
					tagNodeIntensity[bMassMin] += ionInt[i];
				}
/*Fill in the space.*/
				if((bMassMax - bMassMin) > 1)
				{
					for(k = (bMassMin + 1); k < bMassMax; k++)
					{
						tagNode[k] = 1;
						tagNodeIntensity[k] += ionInt[i];
					}
				}
			}
			j--;
		}
	}


	
	return;
}

/**********************FreeTagStructs****************************************
*
* Used for freeing memory in a linked list.  Bob DuBose tells me its best to free 
* space in the reverse order
* that the space was malloc'ed.  This routine does that very thing.
*
*/
void FreeTagStructs(struct Sequence *currPtr)
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
/*****************************MaskSequenceNodeWithTags***************************************
*
*	MaskSequenceNodeWithTags ....blah blah blah
*
*/

void MaskSequenceNodeWithTags(SCHAR *sequenceNode, char *tagNode)
{

	INT_4 i, maxTagRegion;
	INT_4 firstNode, currentNode;
	INT_4 charge = gParam.chargeState - 1;
	INT_4 stop, firstZeroNode;
	REAL_4 precursor;
	tsequence *currentTag;

	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	maxTagRegion = gParam.peptideMW  - 
				   (precursor * charge) - (charge - 1) * gElementMass_x100[HYDROGEN]	/*if precursor was a y ion*/
				   + 2 * gElementMass_x100[HYDROGEN];		/*This is the equation for converting y ion to b ion.*/
				   
				   
	if(gParam.maxent3)
	{
		charge = 1;	/*all maxent3 processed fragment ions are +1*/
	}
	
/*	Initialize tagNode to zero's.*/
	for(i = 0; i < gGraphLength; i++)	
	{
		tagNode[i] = 0;
	}
	
/*	Find the N-terminal node.  This will equal the nominal mass of the N-terminal group R-NH-.*/

	firstNode = gParam.modifiedNTerm;

	tagNode[firstNode] = 1;	/*Give the N-terminal node a value of one.*/

/*	For each tag, assign the nodes a value of one.*/
			
	currentTag = TagSeqList->seq;
	while(currentTag < TagSeqList->seq + TagSeqList->numObjects)
	{
		currentNode = firstNode;
		for(i = 0; i < currentTag->peptideLength; i++)
		{
			currentNode = currentNode + currentTag->peptide[i];
			if(currentNode >= gGraphLength)
			{
				printf("MaskSequenceNodeWithTags:  currentNode >= gGraphLength\n");
				exit(1);
			}
			tagNode[currentNode] = 1;
		}
		currentTag++;
	}
	
/*	
	Find the highest index value in tagNode that is non-zero, add 57 (for Gly)
	and from that node on up reassign each node to a value of one.
*/

	i = gGraphLength - 1;
	while(tagNode[i] == 0)
	{
		i--;
	}
	i = i + gMonoMass_x100[G];
	
/*	If the highest index value plus 57 is less than the maxTagRegion (based on assuming
	only y ions and that no y ions are used that are less than the precursor), then
	maxTagRegion is dropped to the highest index plus 57.
*/

	if(maxTagRegion > i)	
	{
		maxTagRegion = i;
	}
		
	
	for(i = maxTagRegion; i < gGraphLength; i++)
	{
		tagNode[i] = 1;
	}
	
/*	If the N-terminal position contains more than two aa's, then give a value of one to 
	these N-terminal positions in the array.*/
	currentTag = TagSeqList->seq;
	currentNode = firstNode;	/*firstNode is the N-terminal group*/
	while(currentTag < TagSeqList->seq + TagSeqList->numObjects)
	{
		if(currentTag->peptide[0] > currentNode)
		{
			currentNode = currentTag->peptide[0] + firstNode;
		}
		currentTag++;
	}
	if(currentNode >= gGraphLength)
	{
		printf("MaskSequenceNodeWithTags:  currentNode >= gGraphLength\n");
		exit(1);
	}
	for(i = 0; i <= currentNode; i++)
	{
		tagNode[i] = 1;
	}
	
	firstZeroNode = currentNode + 1;
	
/*	For tags that contain W, N, Q, R, or W, the program inserts a value of one at those
	positions corresponding to two amino acid combinations.  This allows data in regions
	outside of the tag region to contribute to the sequencing if there is a gap in the 
	tag region.
*/
	currentTag = TagSeqList->seq;
	while(currentTag < TagSeqList->seq + TagSeqList->numObjects)
	{
		currentNode = firstNode;
		for(i = 0; i < currentTag->peptideLength; i++)
		{
			currentNode += currentTag->peptide[i];
			if(currentNode >= gGraphLength)
			{
				printf("MaskSequenceNodeWithTags:  currentNode >= gGraphLength\n");
				exit(1);
			}
			if(currentTag->peptide[i] == gGapList[N])	/*Asn*/
			{
				if(gGapList[G] != 0)
				{
					tagNode[currentNode - gGapList[G]] = 1;
				}
			}
			if(currentTag->peptide[i] == gGapList[K] 
			   || currentTag->peptide[i] == gGapList[Q])	/*Lys or Gln*/
			{
				if(gGapList[A] != 0)
				{
					tagNode[currentNode - gGapList[A]] = 1;
				}
				if(gGapList[G] != 0)
				{
					tagNode[currentNode - gGapList[G]] = 1;
				}
			}
			if(currentTag->peptide[i] == gGapList[R])	/*Arg*/
			{
				if(gGapList[G] != 0)
				{
					tagNode[currentNode - gGapList[G]] = 1;
				}
				if(gGapList[V] != 0)
				{
					tagNode[currentNode - gGapList[V]] = 1;
				}
			}
			if(currentTag->peptide[i] == gNomMass[W])	/*Trp*/
			{
				if(gGapList[A] != 0)
				{
					tagNode[currentNode - gGapList[A]] = 1;
				}
				if(gGapList[D] != 0)
				{
					tagNode[currentNode - gGapList[D]] = 1;
				}
				if(gGapList[E] != 0)
				{
					tagNode[currentNode - gGapList[E]] = 1;
				}
				if(gGapList[G] != 0)
				{
					tagNode[currentNode - gGapList[G]] = 1;
				}
				if(gGapList[S] != 0)
				{
					tagNode[currentNode - gGapList[S]] = 1;
				}
				if(gGapList[V] != 0)
				{
					tagNode[currentNode - gGapList[V]] = 1;
				}
			}
		}
	
		currentTag++;
	}
	
/*	Assign value of one to gParam.fragmentErr region around each node.*/
	
	i = firstZeroNode;
	while(i < maxTagRegion && i < gGraphLength)
	{
		if(i != firstNode && tagNode[i] != 0)
		{
			i = i - gParam.fragmentErr;
			stop = i + gParam.fragmentErr * 2;
			while(i <= stop && i < gGraphLength)
			{
				tagNode[i] = 1;
				i++;
			}
		}
		i++;
	}
	for(i = firstZeroNode; i <= firstZeroNode + gParam.fragmentErr; i++)
	{
		tagNode[i] = 1;
	}

/*
*	The composite tagNode is multiplied against the array sequenceNode, so that 
*   only sequences identified as possible tags are allowed.
*/

	for(i = 0; i < gGraphLength; i++)
	{
		sequenceNode[i] = tagNode[i] * sequenceNode[i];
	}
		
	return;
}

/*************************************TagExtensions********************************************
*
*	This function adds single amino acid extensions onto the subsequences found in the linked
*	list starting w/ subsequencePtr.  As a result, a new linked list starting with 
*	newSubsequencePtr is generated.  When all of the extensions have been made and the list
*	starting with newSubsequencePtr is complete, then the list starting with subsequencePtr
*	is free'ed, and newSubsequencePtr is returned.  
*	The scores for the subsequences are derived from the array
*	tagNodeIntensity, and are a summation of ion intensities.  The
*	array gGapList and gGapListIndex describe the possible extensions that are allowed.
*	The
*	finalSeqNum is the upper limit on the number of completed sequences that will be stored
*	in the list of Sequence structs in TagSeqList.  The topSeqNum
*	is the maximum number of subsequences that are in the linked list of sequences (both
*	the list starting with subsequencePtr and newSubsequencePtr).
*/

INT_4 TagExtensions(char *tagNode, INT_4 topSeqNum, INT_4 *tagNodeIntensity,
				   INT_4 totalIonIntensity)
{

	tsequence *currPtr;
	tsequence subsequenceToAdd;
	char test;
	INT_4 i, j, testValue;
	INT_4 extNum;
	INT_4 score;
	INT_4 plusPro[AMINO_ACID_NUMBER];
	INT_4 flagToStop = 0;	/*if any subsequence can be extended, then this becomes 1*/
	INT_4 *peptide, peptideLength, nodeValue, gapNum;
	REAL_4 finalScore;
	REAL_4 scoreAdjuster;
	tsequenceList *NewTagSubseqList;
	
	struct extension clearExtension;
	struct extension *extensionList;
	
	clearExtension.gapSize      = 0;
	clearExtension.mass		    = 0;
	clearExtension.singleAAFLAG = 0;
	clearExtension.score        = 0;
	
	extensionList = (extension *) calloc(MAX_GAPLIST, sizeof(struct extension));
	if(extensionList == NULL)
	{
		printf("TagExtensions:  Out of memory.");
		exit(1);
	}
	
	NewTagSubseqList = (tsequenceList *) CreateNewList( sizeof(tsequence), 50, 10 );
	if (!NewTagSubseqList) 
	{
		free(extensionList);
		return(flagToStop);
	}
	
	peptide = (int *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4 ));
	if(peptide == NULL)
	{
		printf("TagExtensions:  Out of memory.");
		exit(1);
	}

/*	Calculate values for two aa extensions that contain proline.*/
	for(i = 0; i < gAminoAcidNumber; i++)
	{
		if(gGapList[i] == 0)
		{
			plusPro[i] = 0;
		}
		else
		{
			plusPro[i] = gGapList[P] + gGapList[i];
		}
	}
	
	currPtr = TagSubseqList->seq;
	while(currPtr < TagSubseqList->seq + TagSubseqList->numObjects)
	{
		/*Clear the extension list*/
		for(i = 0; i < MAX_GAPLIST; i++)
		{
			extensionList[i] = clearExtension;
		}
		
		extNum = 0;
		for(i = 0; i < gAminoAcidNumber; i++)	/*Find the one amino acid extensions.*/
		{
			if(gGapList[i] != 0)
			{
				testValue = currPtr->nodeValue + gGapList[i];
				if(tagNode[testValue] != 0)
				{
					extensionList[extNum].mass         = gGapList[i];
					extensionList[extNum].gapSize      = 0;
					extensionList[extNum].score        = tagNodeIntensity[testValue];
					extensionList[extNum].singleAAFLAG = 1;
					extNum++;
					if(extNum >= MAX_GAPLIST)
					{
						printf("TagExtensions:  extNum >= MAX_GAPLIST\n");
						exit(1);
					}
				}
			}
		}
/*	If no extensions were found, then look for two aa extensions that contain proline.*/
		if(extNum == 0)
		{
			for(i = 0; i < gAminoAcidNumber; i++)	/*Find the one amino acid extensions.*/
			{
				if(plusPro[i] != 0)
				{
					testValue = currPtr->nodeValue + plusPro[i];
					if(tagNode[testValue] != 0)
					{
						extensionList[extNum].mass         = plusPro[i];
						extensionList[extNum].gapSize      = 0;
						extensionList[extNum].score        = tagNodeIntensity[testValue];
						extensionList[extNum].singleAAFLAG = 1;
						extNum++;
						if(extNum >= MAX_GAPLIST)
						{
							printf("TagExtensions:  extNum >= MAX_GAPLIST\n");
							exit(1);
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
		
		if(extNum > 0)
		{
			extensionList = SortExtension(extensionList);
			
			flagToStop = 1;	/*keep extending the subsequences*/
			
		
					
/*	Store this information in the linked list of Sequence structs.  The values placed in the
*	arrays extensionList are used to determine
*	values for the variables peptide[], score, peptideLength, gapNum, and nodeValue.  These
*	variables are passed to a few functions that are used to set up the linked list of structs
*	of type Sequence, which contain the next set of subsequences (newSubsequencePtr*/

			for(j = 0; j < currPtr->peptideLength; j++)	/*Set up the peptide field so that it
														includes the previous values contained
														in the peptide field of currPtr (the
														subsequence currently under investigation.*/
			{
				if(j >= MAX_PEPTIDE_LENGTH)
				{
					printf("TagExtensions:  j >= MAX_PEPTIDE_LENGTH\n");
					exit(1);
				}
				peptide[j] = currPtr->peptide[j];
			}
			
			i = 0;
			while(extensionList[i].mass > 0 && i < MAX_GAPLIST)
			{
				test = TRUE;	/*Becomes FALSE if the data is stored as a final sequence.*/
				score = currPtr->score + extensionList[i].score;
				peptideLength = currPtr->peptideLength + 1;
				peptide[peptideLength - 1] = extensionList[i].mass;
				gapNum = currPtr->gapNum + extensionList[i].gapSize;
				nodeValue = currPtr->nodeValue + extensionList[i].mass;
				if(peptideLength >= MAX_PEPTIDE_LENGTH)
				{
					printf("TagExtensions:  peptideLength >= MAX_PEPTIDE_LENGTH\n");
					exit(1);
				}
				
				/*If there are not too many subsequences stored, then do this.*/
				if(NewTagSubseqList->numObjects < topSeqNum)	
				{
					subsequenceToAdd.score = score;
					for(j = 0; j < peptideLength; j++) 
					{
						subsequenceToAdd.peptide[j] = peptide[j];
					}
					subsequenceToAdd.peptideLength = peptideLength;
					subsequenceToAdd.gapNum = gapNum;
					subsequenceToAdd.nodeValue = nodeValue;
					if(!AddToList(&subsequenceToAdd, NewTagSubseqList)) 
					{
						DisposeList(NewTagSubseqList); 
						free(extensionList);
						free(peptide);
						return(flagToStop);
					}
				}
				else	/*If I have the max allowed number of subsequence, 
									then do this.*/
				{
					qsort(NewTagSubseqList->seq,
			              (size_t)NewTagSubseqList->numObjects,
	   	                  (size_t)sizeof(tsequence),SequenceScoreDescendSortFunc);
					
					if(score > NewTagSubseqList->seq[NewTagSubseqList->numObjects - 1].score)
					{   /* Replace the worst scorer with this sequence. */
						for(j = 0; j < peptideLength; j++)
						{
							NewTagSubseqList->seq[NewTagSubseqList->numObjects -1].peptide[j] = peptide[j];
						}
						NewTagSubseqList->seq[NewTagSubseqList->numObjects -1].peptideLength = peptideLength;
						NewTagSubseqList->seq[NewTagSubseqList->numObjects -1].score = score;
						NewTagSubseqList->seq[NewTagSubseqList->numObjects -1].nodeValue = nodeValue;
						NewTagSubseqList->seq[NewTagSubseqList->numObjects -1].gapNum = gapNum;
					}
				}
				i++;
			}
		}
		if(extNum == 0 && currPtr->peptideLength > 1)
		{
			if(totalIonIntensity == 0)
			{
				printf("TagExtensions:  totalIonIntensity == 0\n");
				exit(1);
			}
			finalScore = (REAL_4)currPtr->score / (REAL_4)totalIonIntensity;
			finalScore = finalScore * 100;
			if((finalScore >= TAG_CUTOFF && 
				(gParam.fragmentPattern == 'Q' || gParam.fragmentPattern == 'T')) ||
				(finalScore >= (TAG_CUTOFF * 0.5) && gParam.fragmentPattern == 'L'))
			{
				for(j = 0; j < currPtr->peptideLength; j++)
				{
					peptide[j] = currPtr->peptide[j];
				}
				score = finalScore;
				peptideLength = currPtr->peptideLength;
				gapNum = currPtr->gapNum;
				nodeValue = currPtr->nodeValue;
				
				/*Here's the expected average peptide length.*/
				scoreAdjuster = nodeValue;				
				scoreAdjuster = scoreAdjuster / gAvResidueMass;
				/*Here's the ratio between the average
				  expected and actual length.*/
				if(peptideLength == 0)
				{
					printf("TagExtensions:  peptideLength == 0\n");
					exit(1);
				}
				scoreAdjuster = scoreAdjuster / peptideLength;	

				if(TagSeqList->numObjects < gParam.finalSeqNum)
				{
					tsequence TagToAdd;
					
					for(j = 0; j < peptideLength; j++)
					{
						TagToAdd.peptide[j] = peptide[j];
					}
					TagToAdd.peptideLength = peptideLength;
					TagToAdd.score = score * scoreAdjuster;
					TagToAdd.nodeValue = nodeValue;
					TagToAdd.gapNum = gapNum;
					if(!AddToList(&TagToAdd, TagSeqList)) 
					{
						DisposeList(NewTagSubseqList);
						free(extensionList);
						free(peptide);
						return(flagToStop);
					}
				}
				else
				{
					qsort(TagSeqList->seq,
			              (size_t)TagSeqList->numObjects,
	   	                  (size_t)sizeof(tsequence),SequenceScoreDescendSortFunc);
					
					if(score > TagSeqList->seq[TagSeqList->numObjects - 1].score)
					{   /* Replace the worst scorer with this sequence. */
						for(j = 0; j < peptideLength; j++)
						{
							TagSeqList->seq[TagSeqList->numObjects -1].peptide[j] = peptide[j];
						}
						TagSeqList->seq[TagSeqList->numObjects -1].peptideLength = peptideLength;
						TagSeqList->seq[TagSeqList->numObjects -1].score = score;
						TagSeqList->seq[TagSeqList->numObjects -1].nodeValue = nodeValue;
						TagSeqList->seq[TagSeqList->numObjects -1].gapNum = gapNum;
					}
				}
			}
		}
		currPtr++;	/*Go to the subsequence.*/
	}

	DisposeList(TagSubseqList); /* Throw away the old list of subsequences */
	TagSubseqList = NewTagSubseqList; /* Make the new subseqs the list of subsequences */
	
	free(extensionList);
	free(peptide);

	return(flagToStop);
}



/*********************************	FindTagYIons***************************************************
*
*	This function assumes that the CID ions are all of type y.  The nominal mass values are 
*	determined and the corresponding positions in the array sequenceNodeC are assigned the 
*	additional value of gWeightedIonValues.y.  
*/

void FindTagYIons(char *tagNode, INT_4 charge, INT_4 *tagNodeIntensity)
{
	tMSData *tagMassPtr;
	INT_4 yMassMin, yMassMax, j;
	REAL_4 yMass, peptideMass;
	REAL_8 aToMFactor;
	
	
	tagMassPtr = TagMassList->mass;
	while(tagMassPtr < TagMassList->mass + TagMassList->numObjects)
	{
		yMass = tagMassPtr->mOverZ;
		yMass = (yMass * charge) - ((charge - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/

		if(yMass >= gParam.monoToAv)	/*Fully apply the average to mono mass factor*/
		{
			aToMFactor = 0;
		}
		else
		{
			if(yMass > (gParam.monoToAv - gAvMonoTransition))	/*In the transition region*/
			{
				aToMFactor = (gParam.monoToAv - yMass) / gAvMonoTransition;
			}
			else	/*Don't apply this factor*/
			{
				aToMFactor = 1;
			}
		}			
		aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
		if(yMass > (gParam.monoToAv - gAvMonoTransition))	/*Convert the y ion to monoisotopic.*/
		{
			yMass = yMass * aToMFactor;
		}
				
		peptideMass = gParam.peptideMW;	
		if(peptideMass >= gParam.monoToAv)	/*Fully apply the average to mono mass factor*/
		{
			aToMFactor = 0;
		}
		else
		{
			if(peptideMass > (gParam.monoToAv - gAvMonoTransition))	/*In the transition region*/
			{
				aToMFactor = (gParam.monoToAv - peptideMass) / gAvMonoTransition;
			}
			else	/*Don't apply this factor*/
			{
				aToMFactor = 1;
			}
		}			
		aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
		if(peptideMass > (gParam.monoToAv - gAvMonoTransition))	/*Convert the y ion to monoisotopic.*/
		{
			peptideMass = peptideMass * aToMFactor;
		}

/*Peptide mass and ymass are converted to monoisotopic before converting to the b ion value*/
		yMass = peptideMass - yMass + (2 * gElementMass_x100[HYDROGEN]);	/*Convert to b ion.*/
							
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node, and convert to the b ion mass.*/
		yMass = yMass + gParam.fragmentErr;
		yMassMax = yMass;		/*Truncate w/o rounding up.*/
			
/*Next I'll find the lowest mass node.*/
		yMass = yMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
		yMassMin = yMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
		if(yMassMax >= gGraphLength)
		{
			printf("FindTagYIons:  yMassMax >= gGraphLength\n");
			exit(1);
		}
		for(j = yMassMin; j <= yMassMax; j++)
		{
			tagNode[j] = 1;
			tagNodeIntensity[j] = tagMassPtr->intensity;
		}
		tagMassPtr++;
	}
	return;
}


/*******************************TagNodeInit************************************************
*
*	This function initializes the array tagNode by first assigning the value of zero
*	to each element in the array (from 0 -> 3999).  Next an N-terminal value is assigned to 
*	position number 1, and then the possible C-terminal nodes are found and given an N-terminal
*	value.  There may be more than one C-terminal node, depending on the mass and the error.
*/

void TagNodeInit(char *tagNode, INT_4 *tagNodeIntensity)
{
	INT_4 i, firstNode;
	REAL_4 lastNode;
	REAL_8 aToMFactor;
	INT_4 lastNodeHigh, lastNodeLow;

	for(i = 0; i < gGraphLength; i++)	/*Initialize tagNode to zero's.*/
	{
		tagNode[i] = 0;
		tagNodeIntensity[i] = 0;
	}
	
/*	Find the N-terminal node.  This will equal the nominal mass of the N-terminal group R-NH-.*/

	firstNode = gParam.modifiedNTerm + 0.5;
	tagNode[firstNode] = 1;	/*Give the N-terminal node a value of one.*/
	
/*Figure out what the C-terminal node(s) are.*/

	lastNode = gParam.peptideMW - gParam.modifiedCTerm;
	
/*Alter the values so that they are closer to the expected nominal masses.*/
	if(lastNode > gParam.monoToAv)
	{
		aToMFactor = 0;
	}
	else
	{	
		if(lastNode > (gParam.monoToAv - gAvMonoTransition))
		{
			aToMFactor = (gParam.monoToAv - lastNode) / gAvMonoTransition;
		}
		else
		{
			aToMFactor = 1;
		}
	}
	aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
	if(lastNode > gParam.monoToAv - gAvMonoTransition)
	{
		lastNode = lastNode * aToMFactor;	/*Convert to monoisotopic*/
	}
	
/*The two extremes for the possible C-terminal nodes are identified first.*/
/*First I'll find the highest mass node.*/
	lastNode = lastNode + gParam.peptideErr;
	lastNodeHigh = lastNode;		/*Truncate w/o rounding up.*/

/*Next I'll find the lowest mass node.*/
	lastNode = lastNode - gParam.peptideErr - gParam.peptideErr + 0.5;
	lastNodeLow = lastNode;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
	if(lastNodeHigh >= gGraphLength)
	{
		printf("TagNodeInit:  lastNodeHigh >= gGraphLength\n");
		exit(1);
	}
	for(i = lastNodeLow; i <= lastNodeHigh; i++)
	{
		tagNode[i] = 1;
	}

	return;
}

/*********************************GetAutoTag******************************************
*
*
*
*/

void GetAutoTag(struct MSData *firstMassPtr, SCHAR *sequenceNode)
{
	struct MSData *currPtr;
	struct Sequence *firstTagPtr;
	INT_4 threshold, charge, *tagNodeIntensity, totalIntensity;
	INT_4 *peptide, count;
	char *tagNode;
	REAL_4 precursor, maxMOverZ;

	TagMassList = (tMSDataList *) CreateNewList( sizeof(tMSData), 10, 10 );
	if (!TagMassList) return;

	TagSubseqList = (tsequenceList *) CreateNewList( sizeof(tsequence), 50, 10 );
	if (!TagSubseqList) 
	{
		DisposeList(TagMassList);
		return;
	}
	
	TagSeqList = (tsequenceList *) CreateNewList( sizeof(tsequence), 10, 10 );
	if (!TagSeqList) 
	{
		DisposeList(TagMassList);
		DisposeList(TagSubseqList);
		return;
	}
	
	
	peptide = (int  *) malloc(MAX_PEPTIDE_LENGTH * sizeof(INT_4));
	if(peptide == NULL)
	{
		printf("GetAutoTag:  Out of memory.");
		exit(1);
	}
	
	tagNode = (char  *) malloc(gGraphLength * sizeof(char));	/*set aside some space for the tagNode*/
	if(tagNode == NULL)
	{
		printf("GetAutoTag:  Out of memory.");
		exit(1);
	}
	tagNodeIntensity = (int *) malloc(gGraphLength * sizeof(INT_4));
	if(tagNodeIntensity == NULL)
	{
		printf("GetAutoTag:  Out of memory.");
		exit(1);
	}

	gTagLength = 0;	/*initialize; used for spectrum quality assessment*/
	
	charge = gParam.chargeState - 1;	/*the only charge considered for the y fragment ions*/
	if(gParam.maxent3)
	{
		charge = 1;	/*all maxent3 processed data is +1*/
	}
	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) / gParam.chargeState;
	maxMOverZ = (((precursor * gParam.chargeState) - gElementMass_x100[HYDROGEN]) / charge)
				 - gGapList[G] + gParam.fragmentErr;
				 
	firstTagPtr = NULL;
	totalIntensity = 0;

	
	if(charge == 0)
	{
		free(tagNode);
		free(tagNodeIntensity);
		free(peptide);
		DisposeList(TagMassList);
		DisposeList(TagSeqList);
		return;	/* Don't try to find tags for singly charged precursors */
	}

/*	
	Find ions that are greater than the precursor but less than the maximum m/z value for
	the given value of charge.  The number 'charge' is the only charge state for the fragment
	ions that is considered here; for +2 precursors I look only for singly charged y ions
	and for +3 precursors I look for +2 y fragment ions.  Ions also must be of intensity 
	greater than 10% of the most abundant ion in this region of the spectrum.
*/

	currPtr = firstMassPtr;	/*find an intensity threshold*/
	count = 0;
	threshold = 0;
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > (precursor + (gParam.fragmentErr * 4))
			&& currPtr->mOverZ < maxMOverZ)
		{
			threshold = threshold + (0.20 * currPtr->intensity);
			count++;
		}
		currPtr = currPtr->next;
	}
	if(count == 0) 
	{
		free(tagNode);
		free(tagNodeIntensity);
		free(peptide);
		DisposeList(TagMassList);
		DisposeList(TagSeqList);
		return;
	}
	threshold = threshold / count;	/*threshold is 20% of the average*/


	/*currPtr = firstMassPtr;*/	/*find 0.1 x most intense ion*/
	/*threshold = 0;
	while(currPtr != NULL)
	{
		if(currPtr->intensity > threshold && currPtr->mOverZ > (precursor + (gParam.peakWidth * 4))
			&& currPtr->mOverZ < maxMOverZ)
		{
			threshold = currPtr->intensity;
		}
		currPtr = currPtr->next;
	}
	threshold = threshold * 0.1;*/	/*Set intensity threshold.*/
	
	/* Find the potential tag ions */
	currPtr = firstMassPtr;	
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > (precursor + (4 * gParam.fragmentErr))
			&& currPtr->mOverZ < maxMOverZ && currPtr != NULL)
		{
			if(currPtr->intensity > threshold)
			{
				tMSData massToAdd;
				
				massToAdd.mOverZ = currPtr->mOverZ;
				massToAdd.intensity = currPtr->intensity;
				massToAdd.normIntensity = currPtr->intensity;
				if(!AddToList(&massToAdd, TagMassList)) 
				{
					free(tagNode);
					free(tagNodeIntensity);
					free(peptide);
					DisposeList(TagMassList);
					DisposeList(TagSeqList);
					return;
				}

				totalIntensity = totalIntensity + currPtr->intensity;	/*set to zero at the top*/
			}
		}
		currPtr = currPtr->next;
	}
	
/*	
*	Filter through the tag masses, so that no ions are closer together than 57/charge.
*	This should not be used for LCQ data, since high mass b ions are mingled with the 
*	high mass y ions.  For LCQ data, I remove neutral loss ions.
*/
	if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q')
	{
		FilterTagMasses(charge);
	}

	if(gParam.fragmentPattern == 'L')
	{
		RemoveNeutralLosses(charge);
	}

/*	
*	Initialize the array tagNode so that all values
*	are zero, except for the N and C terminal nodes, which are assigned a value of one.  
*/
	TagNodeInit(tagNode, tagNodeIntensity);
	
/*
*	Assume all ions in TagMassList are y ions of charge 'charge', and convert these to positions
*	in the tagNode array, where the indexing corresponds to nominal b ion masses.
*/
	FindTagYIons(tagNode, charge, tagNodeIntensity);
	
/*
*	To help reduce the number of sequence tags produced in LCQ data, the b ion nodes are 
*	crossed with the y ion nodes.  That is, only ions that have a corresponding b and y ions
*	are utilized.
*/

	/*if(gParam.fragmentPattern == 'L')
	{
		FindTagBIons(tagNode, charge, tagNodeIntensity);
	}*/
		
/*	
*	Now that the node scores have been finalized (in the array tagNode), 
*	it is time to start building up subsequences from the
*	N-terminus.  I connect the nodes that are spaced one or two amino acid residues apart.
*	The function 
*	TagMaker returns a pointer to a struct of type Sequence (firstTagPtr, 
*	which is the first element in
*	a linked list of completed sequences plus the associated score.  If there were
*	no completed sequences, then the function returns a NULL value.
*/

	/*firstTagPtr = TagMaker(tagNode, tagNodeIntensity, totalIntensity);*/

/*	
	If no tags are found, try sticking in the b2 ion and checking again.
*/
	
	/*if(firstTagPtr == NULL)
	{
		FindTagB2Ions(firstMassPtr, &totalIntensity, tagNode, tagNodeIntensity);
		firstTagPtr = TagMaker(tagNode, tagNodeIntensity, tagMassPtr);
	}*/
	
/*	If still no tags are found, then try searching tagNode from the top down.*/
	if(TagSeqList->numObjects == 0)
	{
		AlternateTagMaker(tagNode, tagNodeIntensity, totalIntensity);
	}
	
/*	
*	If there are any tags found, then these are used to make a composite tagNode, where
*	tagNode array is zeroed out, and reassigned values of one for positions corresponding
*	to the sequences in the tags.  Regions outside of the tag region are assigned values of
*	one.  The composite tagNode is multiplied against the array sequenceNode, so that 
*   only sequences identified as possible tags are allowed when generating subsequences.
*/

	if(TagSeqList->numObjects)
	{
		MaskSequenceNodeWithTags(sequenceNode, tagNode);
	}
	
		
/*	Announce to the world that I'm done w/ the auto-tag.*/
	if(gParam.fMonitor && gCorrectMass)
	{
		if(TagSeqList->numObjects > 1)
		{
			printf("Auto-tag found %ld tags:\n", TagSeqList->numObjects);
		}
		else
		{
			if(TagSeqList->numObjects == 1)
			{
				printf("Auto-tag found one tag:\n");
			}
			else
			{
				printf("Auto-tag found no tags:\n");
			}
		}
		if(TagSeqList->numObjects > 0)
		{   /* Print out the tags */
			tsequence *currentTag;
			
			currentTag = TagSeqList->seq;
			while(currentTag < TagSeqList->seq + TagSeqList->numObjects)
			{
				char *peptideString;
				
				if(currentTag->peptideLength > gTagLength)
				{
					gTagLength = currentTag->peptideLength;
				}
				
				peptideString = ComposePeptideString(currentTag->peptide, 
				                                     currentTag->peptideLength);
				if(peptideString) 
				{
					printf("%s\n", peptideString);
					free(peptideString);
				}

				currentTag++;
			}
		}
	}

/*	Free memory allocations specific to this function.*/
	free(tagNode);
	free(tagNodeIntensity);
	free(peptide);
	DisposeList(TagMassList);
	DisposeList(TagSeqList);
	return;
}
/* -------------------------------------------------------------------------------------
//  SEQUENCE SCORE DESCEND SORT FUNC -- Sort the sequences by their score (High to Low).
*/
int SequenceScoreDescendSortFunc(const void *n1, const void *n2) {
	tsequence *n3, *n4;
	
	n3 = (tsequence *)n1;
	n4 = (tsequence *)n2;
	
	return (int)(n3->score < n4->score)? 1:-1;
}
/* -------------------------------------------------------------------------------------
*  COMPOSE PEPTIDE STRING -- Convert the peptide (stored as masses) into a string,
*                            replacing single amino acids with their one letter equivalent
*                            and masses of ambiguous multiple amino acids with the mass
*                            in brackets.
*                            ex:  peptide = 201, 147, 97  => string = "[201]FP"
*/

char *ComposePeptideString(INT_4 *peptide, INT_4 peptideLength)
{

	INT_4 i, j;
	REAL_4 error = 0.4 * gMultiplier;
	char  test;
	char *string;
	char *p;
	
	
	string = (char *) malloc(50);
	if(!string)
	{
		return NULL;
	}
	p = string;
	
	for(i = 0; i < peptideLength; i++)
	{
		test = TRUE;
		for(j = 0; j < gAminoAcidNumber; j++)
		{
			if(peptide[i] <= gGapList[j] + error && peptide[i] >= gGapList[j] - error)
			{
				p+= sprintf(p, "%c", gSingAA[j]);
				test = FALSE;
				break;
			}
		}
		if(test) /* More than a single AA */
		{
			if(peptide[i] < 1000)
			{
				p+= sprintf(p, "[%3d]", peptide[i]);
			}
			else
			{
				p+= sprintf(p, "[%3d]", peptide[i] / gMultiplier);
			}
		}
	}
	p+= sprintf(p, "\0");	/* NULL terminate the string */
	
	return string;
}	