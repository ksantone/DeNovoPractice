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
	
	LutefiskMakeGraph is a file containing the function MakeSequenceGraph plus its associated
	functions.  It was written to be used as part of a program called "LutefiskXP", which
	is used to aid in the interpretation of CID data of peptides.  The general aim of this file
	(and the function MakeSequenceGraph) is to convert the CID data into a graph of integer
	values corresponding to the nominal masses of singly charged b ions.  There are three
	ways that are being developed for making this conversion (otherwise known as three templates),
	which depend on the type of data that is under investigation.  Currently only one of these
	templates has been completed.  Who can guess when the others might be done?
*/

#include <stdio.h>
#include <stdlib.h>
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"

/********************************RemoveSillyNodes**********************************************
*
*	RemoveSillyNodes removes all nodes below 239 that cannot be made of any combination of 
*	amino acids.  After writing this function I came to doubt that it would really do much
*	good, since the silly nodes are near the N-terminus and won't connect.  Hence, it probably
*	doesn't matter if there are these silly nodes or not.
*/
void RemoveSillyNodes(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN)
{
	INT_4 i, j, k, firstNode, testMass;
	char test;
	
/*	Find the firstNode, based on the N-terminal modification.*/

	firstNode = gParam.modifiedNTerm + 0.5;
	
/*
*	Zero out any non-zero nodes that are less than firstNode + Gly.
*/
	for(i = 0; i < firstNode + gMonoMass_x100[G]; i++)
	{
		if(i != firstNode)
		{
			sequenceNodeC[i] = 0;
			sequenceNodeN[i] = 0;
		}
	}
	
/*
*	Nodes between 57 and less than 142 (2xAla) (plus firstNode) can only be due to 
*	single amino acids.  142 is the smallest number that can be made up from two amino
*	acids that cannot also be made up from one amino acid.  These get zero'ed out.
*/
	for(i = firstNode + gMonoMass_x100[G]; i < firstNode + gMonoMass_x100[A] * 2; i++)
	{
		if(sequenceNodeC[i] != 0 || sequenceNodeN[i] != 0)
		{
			test = TRUE;
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gGapList[j] != 0)
				{
					testMass = firstNode + gGapList[j];
					if(testMass == i)
					{
						test = FALSE;
					}
				}
			}
			if(test)
			{
				sequenceNodeC[i] = 0;
				sequenceNodeN[i] = 0;
			}
		}
	}

/*
*	Step through each node from firstNode + 142 to firstNode + 239.  142 is the smallest number
*	that can be made from two amino acids (2xAla) and 239 (I think) is the smallest number that
*	cannot be made up from two amino acids, but can be made up from three amino acids.
*/
	for(i = firstNode + 142 * gMultiplier; i < firstNode + 239 * gMultiplier; i++)
	{
		if(sequenceNodeC[i] != 0 || sequenceNodeN[i] != 0)
		{
			test = TRUE;
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gGapList[j] != 0)
				{
					for(k = j; k < gAminoAcidNumber; k++)
					{
						if(gGapList[k] != 0)
						{
							testMass = firstNode + gGapList[j] + gGapList[k];
							if(testMass == i)
							{
								test = FALSE;
							}
						}
					}
				}
			}
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(gGapList[j] != 0)
				{
					testMass = firstNode + gGapList[j];
					if(testMass == i)
					{
						test = FALSE;
					}
				}
			}
			if(test)
			{
				sequenceNodeC[i] = 0;
				sequenceNodeN[i] = 0;
			}
		}
	}

	return;
}



/********************************FindTrypticLCQY17Ions********************************************
*
*	This function assumes that the CID ions are all of type y-17 or y-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeC are 
*	assigned the additional value of gWeightedIonValues.y_minus17or18.  Only y-17 or y-18 ions that have a
*	corresponding y ion are counted. 
*/

void FindTrypticLCQY17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeC)
{
	struct MSData *currPtr;
	INT_4 i, j, testForChar;
	INT_4 y17MassMin, y17MassMax;
	REAL_4 y17Mass, peptideMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  
							      This is used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;	/*maxent3 data is converted to +1 fragments*/
	}
	
	/*the max fragment charge state for maxent3 is one*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			y17Mass = currPtr->mOverZ;
			test = IsThisPossible(y17Mass, i);
			
			if(test)
			{
				y17Mass = (y17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
				
				if(y17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(y17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - y17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(y17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					y17Mass = y17Mass * aToMFactor;
				}
				
				peptideMass = gParam.peptideMW;	
				
				if(peptideMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - peptideMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
				{
					peptideMass = peptideMass * aToMFactor;
				}

				
				y17Mass = peptideMass - y17Mass + (2 * gElementMass_x100[HYDROGEN]);	/*Partially convert to 
																	b ion.*/
																				
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				y17Mass = y17Mass - gAmmonia + gParam.fragmentErr;
				y17MassMax = y17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				y17Mass = y17Mass - gParam.fragmentErr - gParam.fragmentErr - gWater + gAmmonia + 0.5;
				y17MassMin = y17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(y17MassMax >= gGraphLength)
				{
					printf("FindTrypticLCQY17Ions:  y17MassMax >= gGraphLength|n");
					exit(1);
				}
				for(j = y17MassMin; j <= y17MassMax; j++)
				{
					if(sequenceNodeC[j] != 0)
					{
						testForChar = sequenceNodeC[j];
						if(i <= mostLikelyFragCharge)
						{
							testForChar += (INT_4)gWeightedIonValues.y_minus17or18;
						}
						else
						{
							testForChar += (INT_4)(gWeightedIonValues.y_minus17or18 * HIGH_CHARGE_MULT);
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeC[j] = testForChar;
						}
						else
						{
							sequenceNodeC[j] = 63;
						}
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}



/*********************************	FindTrypticLCQYIons***************************************************
*
*	This function assumes that the CID ions are all of type y.  The nominal mass values are 
*	determined and the corresponding positions in the array sequenceNodeC are assigned the 
*	additional value of gWeightedIonValues.y.  
*/

void FindTrypticLCQYIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeC)
{
	struct MSData *currPtr;
	INT_4 yMassMin, yMassMax;
	INT_4 i, j, testForChar;
	REAL_4 yMass, peptideMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	
	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  
							      This is used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			yMass = currPtr->mOverZ;
			test = IsThisPossible(yMass, i);
			
			if(test)
			{
				yMass = (yMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
				
				if(yMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(yMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - yMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(yMass > (gParam.monoToAv - gAvMonoTransition))
				{
					yMass = yMass * aToMFactor;
				}
				
				peptideMass = gParam.peptideMW;	
				if(peptideMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - peptideMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
				{
					peptideMass = peptideMass * aToMFactor;
				}

				yMass = peptideMass - yMass + (2 * gElementMass_x100[HYDROGEN]);	/*Convert to b ion.*/
							
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node, and convert to the b ion mass.*/
				if(yMass <= 372 * gMultiplier)	/*Allow greater range for high mass y+2 ions*/
				{
					yMass = yMass + gParam.fragmentErr * i;
				}
				else
				{
					yMass = yMass + gParam.fragmentErr;
				}
				yMassMax = yMass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				if(yMass <= 372 * gMultiplier)	/*Allow greater range for high mass y+2 ions*/
				{
					yMass = yMass - (i * gParam.fragmentErr) - (i * gParam.fragmentErr) + 0.5;
				}
				else
				{
					yMass = yMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				}
				yMassMin = yMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(yMassMax >= gGraphLength)
				{
					printf("FindTrypticLCQYIons:  yMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = yMassMin; j <= yMassMax; j++)
				{
					testForChar = sequenceNodeC[j];	/*make sure value fits in a char*/
					if(i <= mostLikelyFragCharge || yMass <= 373 * gMultiplier)
					{
						testForChar += (INT_4)gWeightedIonValues.y;
					}
					else
					{
						testForChar += (INT_4)(gWeightedIonValues.y * HIGH_CHARGE_MULT);
					}
					if(testForChar < 127 && testForChar > -127)
					{
						sequenceNodeC[j] = testForChar;
					}
					else
					{
						sequenceNodeC[j] = 63;
					}
				}
			}
		}
		currPtr = currPtr->next;
	}
	
/*	If a sequenceNodeC is assigned a non-zero value that is less than the full gWeightedIonValues.y,
*	then that means that it was because an ion was assumed to be of an unlikely charge state.
*	For example, if a doubly charged precursor had an ion that was assumed to be a doubly
*	charged fragment, but the corresponding singly charged ion was absent, then the value
*	for that node would be less than gWeightedIonValues.y.  I remove these from the list here.
*	
*	For the LCQ, I've noticed that there are often doubly charged y ions from +2 precursors, 
*	where the fragment ions are  due to the loss of the N-terminal one, two, or three amino
*	acids.  So don't zero out the first 373 positions (this would be two tryptophans).
*/


	for(i = 373 * gMultiplier; i < gGraphLength; i++)	
	{
		if(sequenceNodeC[i] < gWeightedIonValues.y)
		{
			sequenceNodeC[i] = 0;
		}
	}

	return;
}



/********************************FindTrypticLCQA17Ions********************************************
*
*	This function assumes that the CID ions are all of type a-17 or a-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.a_minus17or18.  Only a-17 or a-18 ions that have a
*	corresponding a ion are counted. 
*/

void FindTrypticLCQA17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent)
{
	struct MSData *currPtr;
	INT_4 a17MassMin, a17MassMax, i, j, testForChar;
	REAL_4 a17Mass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	if(gParam.chargeState == 1)
	{
		mostLikelyFragCharge = 1;
	}
	else
	{	
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			a17Mass = currPtr->mOverZ;
			test = IsThisPossible(a17Mass, i);
			
			if(test)
			{
				a17Mass = (a17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(a17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(a17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - a17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(a17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					a17Mass = a17Mass * aToMFactor;
				}
			
/*The two extremes for the possible a ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				a17Mass = a17Mass + gWater + gCO + gParam.fragmentErr;
				a17MassMax = a17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				a17Mass = a17Mass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				a17MassMin = a17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(a17MassMax >= gGraphLength)
				{
					printf("FindTrypticLCQA17Ions:  a17MassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = a17MassMin; j <= a17MassMax; j++)
				{
					if(ionPresent[j - gCO] != 0)
					{
						testForChar = sequenceNodeN[j];
						if(i <= mostLikelyFragCharge)
						{
							if(a17Mass < 350 * gMultiplier)
							{
								testForChar += (INT_4)gWeightedIonValues.a_minus17or18;
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_MASS_A_MULT);
							}
						}
						else
						{
							if(a17Mass < 350 * gMultiplier)
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_CHARGE_MULT);
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_MASS_A_MULT * HIGH_CHARGE_MULT);
							}
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeN[j] = testForChar;
						}
						else
						{
							sequenceNodeN[j] = 63;
						}
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}


/********************************FindTrypticLCQAIons********************************************
*
*	This function assumes that the CID ions are all of type a.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.a.  Only those a ions that have a
*	corresponding b ion are counted.  Using nominal masses as the index, the number one is 
*	placed in the array ionPresent at the nominal mass of the a ion.
*/

void FindTrypticLCQAIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent)
{
	struct MSData *currPtr;
	INT_4 aMassMin, aMassMax, i, j, testForChar;
	REAL_4 aMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;

	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  This is
							used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			aMass = currPtr->mOverZ;
			test = IsThisPossible(aMass, i);
			
			if(test)
			{
				aMass = (aMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(aMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(aMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - aMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(aMass > (gParam.monoToAv - gAvMonoTransition))
				{
					aMass = aMass * aToMFactor;
				}
			
/*The two extremes for the possible a ion nodes are identified first.  The integer mass value
of the a ion is the index number of the array ionPresent, which is initially set to zero for
all of its values.  If a b ion is present, then the a ion is counted.  The fact that the a ion
has been counted is noted by changing ionPresent[whatever the a ion mass is] to 1.  This is use
when deciding if there are any a-17 ions to be included in the scores for the nodes.
*/
/*First I'll find the highest mass node.*/
				aMass = aMass + gCO + gParam.fragmentErr;
				aMassMax = aMass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				aMass = aMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				aMassMin = aMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(aMassMax >= gGraphLength)
				{
					printf("FindTrypticLCQAIons: aMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = aMassMin; j <= aMassMax; j++)
				{
					if(sequenceNodeN[j] != 0)
					{
						testForChar = sequenceNodeN[j];	/*make sure value fits in a char*/
						if(i <= mostLikelyFragCharge)
						{
							if(aMass < 350 * gMultiplier)
							{
								testForChar += (INT_4)gWeightedIonValues.a;
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_MASS_A_MULT);
							}
						}
						else
						{
							if(aMass < 350 * gMultiplier)
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_CHARGE_MULT);
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_CHARGE_MULT * HIGH_MASS_A_MULT);
							}
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeN[j] = testForChar;
						}
						else
						{
							sequenceNodeN[j] = 63;
						}
						ionPresent[j - 28 * gMultiplier] = 1;	/*The a ion is found.*/
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}



/********************************FindTrypticLCQB17Ions********************************************
*
*	This function assumes that the CID ions are all of type b-17 or b-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.b_minus17or18.  Only b-17 or b-18 ions that have a
*	corresponding b ion are counted. 
*/

void FindTrypticLCQB17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN)
{
	struct MSData *currPtr;
	INT_4 b17MassMin, b17MassMax, i, j, testForChar;
	REAL_4 b17Mass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
														
	if(gParam.chargeState == 1)
	{
		mostLikelyFragCharge = 1;
	}
	else
	{	
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			b17Mass = currPtr->mOverZ;
			test = IsThisPossible(b17Mass, i);
			
			if(test)
			{
				b17Mass = (b17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(b17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(b17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - b17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(b17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					b17Mass = b17Mass * aToMFactor;
				}
							
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				b17Mass = b17Mass + gWater + gParam.fragmentErr;
				b17MassMax = b17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				b17Mass = b17Mass - gParam.fragmentErr - gParam.fragmentErr - gWater + gAmmonia + 0.5;
				b17MassMin = b17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(b17MassMax >= gGraphLength)
				{
					printf("FindTrypticLCQB17Ions: b17MassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = b17MassMin; j <= b17MassMax; j++)
				{
					testForChar = sequenceNodeN[j];	/*make sure value fits in a char*/
					if(sequenceNodeN[j] != 0)
					{
						if(i <= mostLikelyFragCharge)	/*Charge state of the fragment is ok.*/
						{
							testForChar += (INT_4)gWeightedIonValues.b_minus17or18;
						}
						else	/*Charge state of the fragment is probably too high.*/
						{
							testForChar += (INT_4)(gWeightedIonValues.b_minus17or18 * HIGH_CHARGE_MULT);
						}
					}
					if(testForChar < 127 && testForChar > -127)
					{
						sequenceNodeN[j] = testForChar;
					}
					else
					{
						sequenceNodeN[j] = 63;
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}



/*********************************	FindTrypticLCQBIons***************************************************
*
*	This function assumes that the CID ions are all of type b.  The nominal mass values are 
*	determined and the corresponding positions in the array sequenceNodeN are assigned the 
*	additional value of gWeightedIonValues.b.  Rules are for tryptic peptides fragmented in
*	an ion trap.
*/

void FindTrypticLCQBIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN)
{
	struct MSData *currPtr;
	INT_4 bMassMin, bMassMax, i, j, testForChar;
	REAL_4 bMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;

	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  This is
							used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			bMass = currPtr->mOverZ;
			test = IsThisPossible(bMass, i);

			if(test)
			{
				bMass = (bMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(bMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(bMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - bMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(bMass > (gParam.monoToAv - gAvMonoTransition))
				{
					bMass = bMass * aToMFactor;
				}
			
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				bMass = bMass + gParam.fragmentErr;
				bMassMax = bMass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				bMass = bMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				bMassMin = bMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(bMassMax >= gGraphLength)
				{
					printf("FindTrypticLCQBIons:  bMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = bMassMin; j <= bMassMax; j++)
				{
					testForChar = sequenceNodeN[j];	/*test to make sure value fits a char*/
					
					if(i <= mostLikelyFragCharge)
					{
						if(bMass > 147 * gMultiplier)
						{
							testForChar += (INT_4)gWeightedIonValues.b;
						}
					}
					else
					{
						testForChar += (INT_4)(gWeightedIonValues.b * HIGH_CHARGE_MULT);
					}
					if(testForChar < 127 && testForChar > -127)
					{
						sequenceNodeN[j] = testForChar;
					}
					else
					{
						sequenceNodeN[j] = 63;
					}
				}
			}
		}
		currPtr = currPtr->next;
	}
	
/*	If a sequenceNodeN is assigned a non-zero value that is less than the full gWeightedIonValues.b,
*	then that means that it was because an ion was assumed to be of an unlikely charge state.
*	For example, if a doubly charged precursor had an ion that was assumed to be a doubly
*	charged fragment, but the corresponding singly charged ion was absent, then the value
*	for that node would be less than gWeightedIonValues.y.  I remove these from the list here.
*
*	So far, I've not seen alot of multiply charged b ions in LCQ data of tryptic peptides, where
*	the lower charge states (primarily +1 fragments of +2 precursors) is totally absent.
*/
	for(i = 0; i < gGraphLength; i++)	
	{
		if(sequenceNodeN[i] < gWeightedIonValues.b)
		{
			sequenceNodeN[i] = 0;
		}
	}


	return;
}


/*********************************TrypticLCQTemplate***********************************************
*
*	This function takes the initialized sequenceNode and modifies it according to empirical 
*	rules developed from days of experience interpretting CID spectra of multiply-charged
*	tryptic peptides using the LCQ.  The rules are as follows:
*		Assume each ion is both a b and a y ion.  Look for an a ion if a b ion is present.
*	Look for b-17 ion if a b ion is present.  Look for an a-17 ion if an a ion is present
*	(and the corresponding b ion is also present).  Look for a y-17 ion if a y ion is present.
*	Look for multiply-charged b and y ions if sufficient mass is available.
*	
*/

void TrypticLCQTemplate(struct MSData *firstMassPtr, SCHAR *sequenceNodeC, 
						SCHAR *sequenceNodeN)
{
	char *ionPresent;
	INT_4 i;
	
	ionPresent = (char *) malloc(gGraphLength * sizeof(char ));	/*Will contain C-terminal evidence.*/
	if(ionPresent == NULL)
	{
		printf("TrypticLCQTemplate:  Out of memory");
		exit(1);
	}

	for(i = 0; i < gGraphLength; i++)	/*Init this array to zero.  These arrays are used to 
										keep track of where the a ions are.*/
	{
		ionPresent[i] = 0;
	}
	
	if(gWeightedIonValues.b != 0)
	{
		FindTrypticLCQBIons(firstMassPtr, sequenceNodeN);
	}
	
	if(gWeightedIonValues.b_minus17or18 != 0)
	{
		FindTrypticLCQB17Ions(firstMassPtr, sequenceNodeN);
	}
	
	if(gWeightedIonValues.a != 0)
	{				
		FindTrypticLCQAIons(firstMassPtr, sequenceNodeN, ionPresent);
	}
	
	if(gWeightedIonValues.a_minus17or18	!= 0)
	{				
		FindTrypticLCQA17Ions(firstMassPtr, sequenceNodeN, ionPresent);
	}
	
	if(gWeightedIonValues.y != 0)
	{					
		FindTrypticLCQYIons(firstMassPtr, sequenceNodeC);
	}
	
	if(gWeightedIonValues.y_minus17or18 != 0)
	{			
		FindTrypticLCQY17Ions(firstMassPtr, sequenceNodeC);
	}
	
	free(ionPresent);
	
	return;
}


/****************************RatchetIt************************************
*
*	RatchetIt increases sequence[cycle] by one.  If the value of sequence[cycle] exceeds
*	the aaNum for that cycle position, then sequence[cycle] is reset to zero and cycle is 
*	increased by one and seqeunce[cycle] (old cycle value + 1) is increased by one.  The
*	input is 'aaNum' (the number of amino acids in the current cycle), cycle, the array 'sequence'
*	and 'seqLength' (the maximum length of this sequence).  It returns a value of TRUE until
*	'cycle' exceeds the 'seqLength' (at which point it returns a FALSE char).  By way 
*	of example, the following shows how RatchetIt should handle a list containing three 
*	amino acids per cycle and is two cycles INT_4:
*		0	0
*		1	0
*		2	0
*		0	1
*		1	1
*		2	1
*		0	2
*		1	2
*		2	2	--->		test = 0 (exit the while loop)
*/

char RatchetIt (INT_4 *aaNum, char cycle, char *sequence, INT_4 seqLength)
{
	char test = TRUE;
/*
*	"cycle" is initialized to zero prior to calling RatchetIt for the first time, and is 
*	increased only if the zero index of sequence exceeds the number of amino acids in the 
*	zero-th Edman cycle (which is really the first cycle of Edman degradation.  
*	"sequence[cycle]" is always increased by one, but if it exceeds the number of amino acids 
*	in that "cycle" then it is reset to zero and cycle is increased by one.
*/
	sequence[cycle] += 1;	/*Sometimes, this is all that RatchetIt does.*/
	
	if(sequence[cycle] >= aaNum[cycle] && cycle <= seqLength)	/*Sometimes, it does more.*/
	{
		sequence[cycle] = 0;	/*Reset the array 'sequence' to zero at the current 'cycle'.*/
		cycle += 1;	/*Go to the next highest cycle.*/
		if(cycle > seqLength)	/*If I've gone too far w/ 'cycle', then shut the whole thing down
								by returning a FALSE char value.*/
		{
			test = FALSE;
			return(test);
		}
		test = RatchetIt(aaNum, cycle, sequence, seqLength);	/*Recursive call.*/
	}

	return(test);
}


/*******************************AddEdmanData***************************************************
*
*	This function alters the arrays sequenceNodeC and sequenceNodeN so as to incorporate the
*	additional data provided from Edman sequencing.  All permutations of the Edman data are 
*	determined, and the corresponding nodes (assuming the N-terminus is unmodified) are altered.
*	If its not zero, then a value corresponding to one-half of totalIonVal is added to the
*	existing node value.
*/	

void AddEdmanData(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal)
{
	INT_4 i, j, aaNum[MAX_PEPTIDE_LENGTH], edmanNode, halfTotalIonVal;
	char sequence[MAX_PEPTIDE_LENGTH], cycle;
	
	halfTotalIonVal = totalIonVal * 0.5;	/*This is added to the nodes.*/
	
/*Determine the number of amino acids in each cycle.*/

	for(i = 0; i < gMaxCycleNum; i++)	
	{
		aaNum[i] = 0;
		j = 0;
		while(gEdmanData[i][j] != 0)	/*If there's a zero, then that signals the end of the list.*/
		{
			aaNum[i]++;
			j++;
		}
	}
	
/*	
*	Here's the major loop in this function - it increments through each Edman cycle that was 
*	entered in the file Lutefisk.edman.  What I'm doing here is finding all combinations of 
*	Edman-derived amino acids of varying lengths.  So I start out w/ i = 0, which is all
*	combinations that are one amino acid INT_4 (ie, just the first cycle).  Next I look for i = 1,
*	which is all combinations from the first and second cycle to give a two amino acid INT_4
*	segment.  This continues until I get to the maximum number of Edman cycles available.
*/
	
	for(i = 0; i < gMaxCycleNum; i++)	/*Nodes will contain from one to gMaxCycleNum-1 amino
										acids.*/
	{
/*	Initialize some variables for each time through this loop.*/
		for(j = 0; j < MAX_PEPTIDE_LENGTH; j++)
		{
			sequence[j] = 0;
		}
		sequence[0] = -1;
		cycle = 0;
		
/*	
*	RatchetIt will alter the values found in the array "sequence".  This array corresponds to 
*	the indexing of one of the dimensions in the array gEdmanData[i][sequence[i]].  If cycle 
*	reaches a value of gMaxCycleNum, then it returns a char value of zero; otherwise its one.
*/
		while(RatchetIt(aaNum, cycle, sequence, i))
		{

/*	Calculate the node for this particular Edman-derived sequence.*/

			edmanNode = gElementMass_x100[HYDROGEN];
			for(j = 0; j <= i; j++)
			{
				edmanNode += gEdmanData[j][sequence[j]];
			}
			
/*	If that node is non-zero for 'sequenceNodeC', then add the value 'halfTotalIonVal' to it.*/

			if(sequenceNodeC[edmanNode] != 0)
			{
				sequenceNodeC[edmanNode] += halfTotalIonVal;
			}
			else
			{
				if(i == 0)
				{
					sequenceNodeC[edmanNode] = 1;	/*If its one amino acid INT_4 and zero value,
													then do this.*/
				}
				else
				{
					sequenceNodeC[edmanNode] = 0;	/*If its over one amino acid INT_4 and zero 
													value, then do this.  Currently it does 
													nothing.*/
				}
			}
			if(sequenceNodeN[edmanNode] != 0)	/*Same as above, but for sequenceNodeN.*/
			{
				sequenceNodeN[edmanNode] += halfTotalIonVal;
			}
			else
			{
				if(i == 0)
				{
					sequenceNodeN[edmanNode] = 1;
				}
				else
				{
					sequenceNodeN[edmanNode] = 0;
				}

			}
		}
	}

	return;
}

/***************************************AddCTermResidue****************************************
*
*	This function locates all of the C-terminal Nodes and subtracts the nominal mass of the
*	most likely C-terminal amino acids for a given proteolysis.  If that node has a value of zero,
*	then a value of one is placed there, so that it might be used for subsequence building.
*/

void AddCTermResidue(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN)
{
	INT_4 i = gGraphLength - 1;
	char test = TRUE;

/*	
*	The INT_4 'i' is used to index the arrays 'sequenceNodeC' and 'sequenceNodeN'.  Obviously,
*	it must not be a negative number.  The char 'test' is initialized as TRUE, and becomes
*	FALSE once a non-zero value of sequenceNodeC[i] is encountered (as i is incremented down from
*	its maximum value - GRAPH_LENGTH.  Once 'test' becomes FALSE, then this while loop can 
*	continue only as INT_4 as sequenceNodeC[i] remains non-zero.  Remember that due to mass 
*	measurement errors, there can be several C-terminal nodes, but that they will always be 
*	adjacent and separate from any other nodes by a string of zero-valued nodes.
*/

	while(i > 0 && (test || sequenceNodeC[i] != 0))
	{
		i--;
		if(sequenceNodeC[i] != 0)
		{
			test = FALSE;
			if(gParam.proteolysis == 'T')	/*If its a tryptic cleavage.*/
			{
			/*If LCQ data, give boost to b ions w/ loss of K or R*/
				if(gParam.fragmentPattern == 'L')
				{
					sequenceNodeN[i - gMonoMass_x100[K]] = sequenceNodeN[i - gMonoMass_x100[K]] * 4;
					sequenceNodeN[i - gMonoMass_x100[R]] = sequenceNodeN[i - gMonoMass_x100[R]] * 4;
				}
				else	/*if qtof or triple quad, boost the y2 ions*/
				{
					sequenceNodeC[i - gMonoMass_x100[K]] = sequenceNodeC[i - gMonoMass_x100[K]] * 4;
					sequenceNodeC[i - gMonoMass_x100[R]] = sequenceNodeC[i - gMonoMass_x100[R]] * 4;
				}
				
				if(sequenceNodeC[i - gMonoMass_x100[K]] == 0 && sequenceNodeN[i - gMonoMass_x100[K]] == 0)	/*Lys*/
				{
					sequenceNodeC[i - gMonoMass_x100[K]] = 10;
					sequenceNodeN[i - gMonoMass_x100[K]] = 10;
				}
				if(sequenceNodeC[i - gMonoMass_x100[R]] == 0 && sequenceNodeN[i - gMonoMass_x100[R]] == 0)	/*Arg*/
				{
					sequenceNodeC[i - gMonoMass_x100[R]] = 10;
					sequenceNodeN[i - gMonoMass_x100[R]] = 10;
				}
			}
			if(gParam.proteolysis == 'K')	/*If its a Lys-C cleavage.*/
			{
				if(sequenceNodeC[i - gMonoMass_x100[K]] == 0 && sequenceNodeN[i - gMonoMass_x100[K]] == 0)	/*Lys*/
				{
					sequenceNodeC[i - gMonoMass_x100[K]] = 1;
					sequenceNodeN[i - gMonoMass_x100[K]] = 1;
				}
			}
			if(gParam.proteolysis == 'E')	/*If its a Staph v8 cleavage.*/
			{
				if(sequenceNodeC[i - gMonoMass_x100[E]] == 0 && sequenceNodeN[i - gMonoMass_x100[E]] == 0)	/*Glu*/
				{
					sequenceNodeC[i - gMonoMass_x100[E]] = 1;
					sequenceNodeN[i - gMonoMass_x100[E]] = 1;
				}
				if(sequenceNodeC[i - gMonoMass_x100[D]] == 0 && sequenceNodeN[i - gMonoMass_x100[D]] == 0)	/*Asp*/
				{
					sequenceNodeC[i - gMonoMass_x100[D]] = 1;
					sequenceNodeN[i - gMonoMass_x100[D]] = 1;
				}
			}
		}
	}
	if(gParam.proteolysis == 'D')	/*If its an Asp-N cleavage.*/
	{
		if(sequenceNodeC[gMonoMass_x100[D] + gElementMass_x100[HYDROGEN]] == 0 && sequenceNodeN[gMonoMass_x100[D] + gElementMass_x100[HYDROGEN]] == 0)
		{
			sequenceNodeC[gMonoMass_x100[D] + gElementMass_x100[HYDROGEN]] = 1;
			sequenceNodeN[gMonoMass_x100[D] + gElementMass_x100[HYDROGEN]] = 1;
		}
	}
	return;
}



/*****************************AddTag**********************************************************
*
*	This function inputs the sequence tag (a character string of single letter code amino
*	acids), and the unsequenced mass at the N- and C-terminii (both REAL_4s).  In addition,
*	the mono to average mass switch mass, the peptideMW, fragmentErr, peptideErr are input.  
*	Also, the arrays sequenceNodeC and sequenceNodeN are used to figure out which nodes
*	contain any sequence evidence.  There is no output, except that these two arrays
*	(sequenceNodeN and sequenceNodeC) are modified so that the interveneing sequence tag
*	region is cut out and replaced by  "superNodes".  The superNode values are the sum of the
*	intervening node values that correspond to the sequence tag.  The value of peptideMW is altered
*	so that it is reduced by the mass of the sequence tag.  The modified arrays continue on 
*	to be processed in the usual manner, but when a completed sequence is reached, then the
*	tag is re-inserted in the appropriate location.
*/

void AddTag(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN)
{
	INT_4 minSuperNode, maxSuperNode, i, j, k, minCTerm, maxCTerm;
	INT_4 nextNode, oldValue, newValue;
	REAL_4 superNodeMass, cTermMass;
	REAL_8 aToMFactor;
	char test;
	
/*Convert the C-terminal tag mass to monoisotopic mass.*/

	cTermMass = gParam.tagCMass;	
	if(cTermMass >= gParam.monoToAv)
	{
		aToMFactor = 0;
	}
	else
	{
		if(cTermMass > (gParam.monoToAv - gAvMonoTransition))
		{
			aToMFactor = (gParam.monoToAv - cTermMass) / gAvMonoTransition;
		}
		else
		{
			aToMFactor = 1;
		}
	}
	aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
	if(cTermMass > (gParam.monoToAv - gAvMonoTransition))
	{
		cTermMass = cTermMass * aToMFactor;
	}
	
/*	
*	The C-terminal mass 'tagCMass' includes the mass of the residues C-terminal to the 
*	'tagSequence' plus the C-terminal group (CO-R), which is either -OH or -NH2.  So I
*	need to change the C-terminal mass to reflect the type of C-terminus.  This will make 
*	it so that the thus modified C-terminal mass added to the N-terminal mass plus the masses 
*	of the sequence tag will equal the C-terminal node values (which are sort of like b ions).
*/
	cTermMass = cTermMass - gParam.modifiedCTerm;

/*	Find the maximum c terminal mass.*/
	cTermMass = cTermMass + gParam.fragmentErr;
	maxCTerm = cTermMass;
/*	Find the minimum c terminal mass.*/
	cTermMass = cTermMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
	minCTerm = cTermMass;
	
	superNodeMass = gParam.tagNMass;
	if(superNodeMass >= gParam.monoToAv)
	{
		aToMFactor = 0;
	}
	else
	{
		if(superNodeMass > (gParam.monoToAv - gAvMonoTransition))
		{
			aToMFactor = (gParam.monoToAv - superNodeMass) / gAvMonoTransition;
		}
		else
		{
			aToMFactor = 1;
		}
	}
	aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
	if(superNodeMass > (gParam.monoToAv - gAvMonoTransition))
	{
		superNodeMass = superNodeMass * aToMFactor;
	}
	
/*	Find the maximum super node mass.*/
	superNodeMass = superNodeMass + gParam.fragmentErr;
	maxSuperNode = superNodeMass;
/*	Find the minimum super node mass.*/
	superNodeMass = superNodeMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
	minSuperNode = superNodeMass;
	if(minSuperNode < gElementMass_x100[HYDROGEN])	/*If tag includes N-terminal amino acids,
													but mass is off a bit*/
	{
		minSuperNode = gElementMass_x100[HYDROGEN];
	}
	
	for(i = minSuperNode; i <= maxSuperNode; i++)
	{
		nextNode = i;	/*Start with the node for the N-terminal mass.*/
		j = 0;

/*	
*	Calculate the node (or mass) of the sequence tag plus the N-terminal mass.  Set test to 
*	FALSE, and add to the superNode positions of sequenceNodeN and sequenceNodeC the values
*	found in the nodes that make up the sequence tag.
*/

		while(gParam.tagSequence[j] != 0)
		{
			test = TRUE;
			for(k = 0; k < gAminoAcidNumber; k++)
			{
				if(gGapList[k] != 0)
				{
					if(gSingAA[k] == gParam.tagSequence[j])
					{
						test = FALSE;
						nextNode = nextNode + gGapList[k];
					}
				}
			}
				
			if(test)	/*TRUE only if a sequence tag entry was not a valid amino acid.*/
			{
				printf("There is something wrong with the sequence tag.");
				exit(1);
			}
			j++;
		}
		test = FALSE;
		for(j = nextNode + minCTerm; j <= nextNode + maxCTerm; j++)
		{
			if(sequenceNodeN[j] != 0)
			{
				test = TRUE;
				break;
			}
		}
		
/*	Change the supernode values to a negative one.*/
		if(test)
		{
			if(i >= gGraphLength)
			{
				printf("AddTag:  i >= gGraphLength\n");
				exit(1);
			}
			sequenceNodeN[i] = -1;
			sequenceNodeC[i] = -1;
		}
	}
	
/*	Now I reassign the node values so that the sequence tag region is cut out.*/
	for(i = gGraphLength - 1; i > 0; i--)	/*Find the next node above the superNode series*/
	{
		if(sequenceNodeN[i] == -1)
		{
			oldValue = i + 1;
			break;
		}
	}
	
	j = 0;	/*Find the node that is the mass of the tagSequence higher than the oldValue*/
	newValue = oldValue;
	while(gParam.tagSequence[j] != 0)
	{
		test = TRUE;
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(gGapList[k] != 0)
			{
				if(gSingAA[k] == gParam.tagSequence[j])
				{
					test = FALSE;
					newValue = newValue + gGapList[k];
				}
			}
		}
		if(test)	/*TRUE only if a sequence tag entry was not a valid amino acid.*/
		{
			printf("There is something wrong with the sequence tag.");
			exit(1);
		}
		j++;
	}

	while(newValue < gGraphLength)
	{
		sequenceNodeN[oldValue] = sequenceNodeN[newValue];
		sequenceNodeC[oldValue] = sequenceNodeC[newValue];
		newValue++;
		oldValue++;
	}
	while(oldValue < gGraphLength)
	{
		sequenceNodeN[oldValue] = 0;
		sequenceNodeC[oldValue] = 0;
		oldValue++;
	}
	
/*	Reset the peptideMW to reflect the removal of the sequence tag mass.*/
	j = 0;
	while(gParam.tagSequence[j] != 0)
	{
		for(k = 0; k < gAminoAcidNumber; k++)
		{
			if(gGapList[k] != 0)
			{
				if(gSingAA[k] == gParam.tagSequence[j])
				{
					gParam.peptideMW = gParam.peptideMW - gGapList[k];
				}
			}
		}
		j++;
	}


	return;
}
/********************************FindTrypticY17Ions********************************************
*
*	This function assumes that the CID ions are all of type y-17 or y-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeC are 
*	assigned the additional value of gWeightedIonValues.y_minus17or18.  Only y-17 or y-18 ions that have a
*	corresponding y ion are counted. 
*/

void FindTrypticY17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeC)
{
	struct MSData *currPtr;
	INT_4 i, j, testForChar;
	INT_4 y17MassMin, y17MassMax;
	REAL_4 y17Mass, peptideMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  
							      This is used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			y17Mass = currPtr->mOverZ;
			test = IsThisPossible(y17Mass, i);
			
			if(test)
			{
				y17Mass = (y17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
				
				if(y17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(y17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - y17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(y17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					y17Mass = y17Mass * aToMFactor;
				}
				
				peptideMass = gParam.peptideMW;	
				
				if(peptideMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - peptideMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
				{
					peptideMass = peptideMass * aToMFactor;
				}

				
				y17Mass = peptideMass - y17Mass + (2 * gElementMass_x100[HYDROGEN]);	/*Partially convert to 
																	b ion.*/
																	
			
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				y17Mass = y17Mass - gAmmonia + gParam.fragmentErr;
				y17MassMax = y17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				y17Mass = y17Mass - gParam.fragmentErr - gParam.fragmentErr - gWater + gAmmonia + 0.5;
				y17MassMin = y17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(y17MassMax >= gGraphLength)
				{
					printf("FindTrypticY17Ions: y17MassMax >= gGraphLength\n"); 
					exit(1);
				}
				for(j = y17MassMin; j <= y17MassMax; j++)
				{
					if(sequenceNodeC[j] != 0)
					{
						testForChar = sequenceNodeC[j];
						if(i <= mostLikelyFragCharge)
						{
							testForChar += (INT_4)gWeightedIonValues.y_minus17or18;
						}
						else
						{
							testForChar += (INT_4)(gWeightedIonValues.y_minus17or18 * HIGH_CHARGE_MULT);
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeC[j] = testForChar;
						}
						else
						{
							sequenceNodeC[j] = 63;
						}
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}


/*********************************	FindTrypticYIons***************************************************
*
*	This function assumes that the CID ions are all of type y.  The nominal mass values are 
*	determined and the corresponding positions in the array sequenceNodeC are assigned the 
*	additional value of gWeightedIonValues.y.  
*/

void FindTrypticYIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeC)
{
	struct MSData *currPtr;
	INT_4 yMassMin, yMassMax;
	INT_4 i, j, testForChar, firstNode, massDiff;
	REAL_4 yMass, peptideMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	
	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  
							      This is used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			yMass = currPtr->mOverZ;
			test = IsThisPossible(yMass, i);
			
			if(test)
			{
				yMass = (yMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
				
				if(yMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(yMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - yMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(yMass > (gParam.monoToAv - gAvMonoTransition))
				{
					yMass = yMass * aToMFactor;
				}
				
				peptideMass = gParam.peptideMW;	
				if(peptideMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - peptideMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(peptideMass > (gParam.monoToAv - gAvMonoTransition))
				{
					peptideMass = peptideMass * aToMFactor;
				}

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
					printf("FindTrypticYIons:  yMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = yMassMin; j <= yMassMax; j++)
				{
					testForChar = sequenceNodeC[j];	/*make sure value fits in a char*/
					if(i <= mostLikelyFragCharge)
					{
						testForChar += (INT_4)gWeightedIonValues.y;
					}
					else
					{
						testForChar += (INT_4)(gWeightedIonValues.y * HIGH_CHARGE_MULT);
					}
					if(testForChar < 127 && testForChar > -127)
					{
						sequenceNodeC[j] = testForChar;
					}
					else
					{
						sequenceNodeC[j] = 63;
					}
				}
			}
		}
		currPtr = currPtr->next;
	}
	
/*	If a sequenceNodeC is assigned a non-zero value that is less than the full gWeightedIonValues.y,
*	then that means that it was because an ion was assumed to be of an unlikely charge state.
*	For example, if a doubly charged precursor had an ion that was assumed to be a doubly
*	charged fragment, but the corresponding singly charged ion was absent, then the value
*	for that node would be less than gWeightedIonValues.y.  I remove these from the list here.
*/
	firstNode = 0;
	for(i = 0; i < gGraphLength; i++)	
	{
		if(firstNode == 0)
		{
			if(sequenceNodeC[i] != 0)
			{
				firstNode = i;	/*find the N-terminus*/
			}
		}
		/*if +2 y ion present indicating N-terminal amino acid, then it will have
		a node value less than the max*/
		if(sequenceNodeC[i] > 0 && sequenceNodeC[i] < gWeightedIonValues.y && i != firstNode)
		{
			massDiff = i - firstNode;
			for(j = 0; j < gAminoAcidNumber; j++)
			{
				if(massDiff == gGapList[j])
				{
					sequenceNodeC[i] = gWeightedIonValues.y;	/*give it higher value so that
																it can be retained below*/
				}
			}
		}
		if(sequenceNodeC[i] < gWeightedIonValues.y)
		{
			sequenceNodeC[i] = 0;
		}
	}

	return;
}

/********************************FindTrypticA17Ions********************************************
*
*	This function assumes that the CID ions are all of type a-17 or a-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.a_minus17or18.  Only a-17 or a-18 ions that have a
*	corresponding a ion are counted. 
*/

void FindTrypticA17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent)
{
	struct MSData *currPtr;
	INT_4 a17MassMin, a17MassMax, i, j, testForChar;
	REAL_4 a17Mass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	if(gParam.chargeState == 1)
	{
		mostLikelyFragCharge = 1;
	}
	else
	{	
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			a17Mass = currPtr->mOverZ;
			test = IsThisPossible(a17Mass, i);
			
			if(test)
			{
				a17Mass = (a17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(a17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(a17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - a17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(a17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					a17Mass = a17Mass * aToMFactor;
				}
			
/*The two extremes for the possible a ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				a17Mass = a17Mass + gWater + gCO + gParam.fragmentErr;
				a17MassMax = a17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				a17Mass = a17Mass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				a17MassMin = a17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(a17MassMax >= gGraphLength)
				{
					printf("FindTrypticA17Ions:  a17MassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = a17MassMin; j <= a17MassMax; j++)
				{
					if(ionPresent[j - 28 * gMultiplier] != 0)
					{
						testForChar = sequenceNodeN[j];
						if(i <= mostLikelyFragCharge)
						{
							if(a17Mass < 350 * gMultiplier)
							{
								testForChar += (INT_4)gWeightedIonValues.a_minus17or18;
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_MASS_A_MULT);
							}
						}
						else
						{
							if(a17Mass < 350 * gMultiplier)
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_CHARGE_MULT);
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a_minus17or18 * HIGH_MASS_A_MULT * HIGH_CHARGE_MULT);
							}
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeN[j] = testForChar;
						}
						else
						{
							sequenceNodeN[j] = 63;
						}
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}

/********************************FindTrypticAIons********************************************
*
*	This function assumes that the CID ions are all of type a.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.a.  Only those a ions that have a
*	corresponding b ion are counted.  Using nominal masses as the index, the number one is 
*	placed in the array ionPresent at the nominal mass of the a ion.
*/

void FindTrypticAIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent)
{
	struct MSData *currPtr;
	INT_4 aMassMin, aMassMax, i, j, testForChar;
	REAL_4 aMass;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;

	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  This is
							used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			aMass = currPtr->mOverZ;
			test = IsThisPossible(aMass, i);
			
			if(test)
			{
				aMass = (aMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(aMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(aMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - aMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(aMass > (gParam.monoToAv - gAvMonoTransition))
				{
					aMass = aMass * aToMFactor;
				}
			
/*The two extremes for the possible a ion nodes are identified first.  The integer mass value
of the a ion is the index number of the array ionPresent, which is initially set to zero for
all of its values.  If a b ion is present, then the a ion is counted.  The fact that the a ion
has been counted is noted by changing ionPresent[whatever the a ion mass is] to 1.  This is use
when deciding if there are any a-17 ions to be included in the scores for the nodes.
*/
/*First I'll find the highest mass node.*/
				aMass = aMass + gCO + gParam.fragmentErr;
				aMassMax = aMass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				aMass = aMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				aMassMin = aMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(aMassMax >= gGraphLength)
				{
					printf("FindTrypticAIons:  aMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = aMassMin; j <= aMassMax; j++)
				{
					if(sequenceNodeN[j] != 0)
					{
						testForChar = sequenceNodeN[j];	/*make sure value fits in a char*/
						if(i <= mostLikelyFragCharge)
						{
							if(aMass < 350 * gMultiplier)
							{
								testForChar += (INT_4)gWeightedIonValues.a;
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_MASS_A_MULT);
							}
						}
						else
						{
							if(aMass < 350 * gMultiplier)
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_CHARGE_MULT);
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.a * HIGH_CHARGE_MULT * HIGH_MASS_A_MULT);
							}
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeN[j] = testForChar;
						}
						else
						{
							sequenceNodeN[j] = 63;
						}
						ionPresent[j - 28 * gMultiplier] = 1;	/*The a ion is found.*/
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}

/********************************FindTrypticB17Ions********************************************
*
*	This function assumes that the CID ions are all of type b-17 or b-18.  The nominal mass 
*	values are determined and the corresponding positions in the array sequenceNodeN are 
*	assigned the additional value of gWeightedIonValues.b_minus17or18.  Only b-17 or b-18 ions that have a
*	corresponding b ion are counted. 
*/

void FindTrypticB17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN)
{
	struct MSData *currPtr;
	INT_4 b17MassMin, b17MassMax, i, j, testForChar;
	REAL_4 b17Mass, precursor;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	precursor = gParam.peptideMW / gParam.chargeState;	/*Not exactly the precursor, but close 
														enough.*/
														
	if(gParam.chargeState == 1)
	{
		mostLikelyFragCharge = 1;
	}
	else
	{	
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			b17Mass = currPtr->mOverZ;
			test = IsThisPossible(b17Mass, i);
			
			if(test)
			{
				b17Mass = (b17Mass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(b17Mass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(b17Mass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - b17Mass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(b17Mass > (gParam.monoToAv - gAvMonoTransition))
				{
					b17Mass = b17Mass * aToMFactor;
				}
			
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				b17Mass = b17Mass + gWater + gParam.fragmentErr;
				b17MassMax = b17Mass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				b17Mass = b17Mass - gParam.fragmentErr - gParam.fragmentErr - gWater + gAmmonia + 0.5;
				b17MassMin = b17Mass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(b17MassMax >= gGraphLength)
				{
					printf("FindTrypticB17Ions:  b17MassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = b17MassMin; j <= b17MassMax; j++)
				{
					testForChar = sequenceNodeN[j];	/*make sure value fits in a char*/
					if(sequenceNodeN[j] != 0)
					{
						if(i <= mostLikelyFragCharge)	/*Charge state of the fragment is ok.*/
						{
							if(b17Mass < precursor)	/*If the m/z is less than the precursor, thats good.*/
							{
								testForChar += (INT_4)gWeightedIonValues.b_minus17or18;
							}
							else	/*Otherwise, there's a penalty.*/
							{
								testForChar += (INT_4)(gWeightedIonValues.b_minus17or18 * HIGH_MASS_B_MULT);
							}
						}
						else	/*Charge state of the fragment is probably too high.*/
						{
							if(b17Mass < precursor)
							{
								testForChar += (INT_4)(gWeightedIonValues.b_minus17or18 * HIGH_CHARGE_MULT);
							}
							else
							{
								testForChar += (INT_4)(gWeightedIonValues.b_minus17or18 * HIGH_CHARGE_MULT * HIGH_MASS_B_MULT);
							}
						}
						if(testForChar < 127 && testForChar > -127)
						{
							sequenceNodeN[j] = testForChar;
						}
						else
						{
							sequenceNodeN[j] = 63;
						}
					}
				}
			}
		}
		currPtr = currPtr->next;
	}

	return;
}

/********************************IsThisStillPossible*******************************************
*
*	This function is called if an ion under consideration is of greater m/z than the precursor.
*	If there is a possible a ion present at 28/charge below the bMass, then a TRUE is returned;
*	otherwise a FALSE is returned.
*
*/
char IsThisStillPossible(REAL_4 bMass, INT_4 currentCharge, struct MSData *firstMassPtr)
{
	char test = FALSE;
	REAL_4 testMass;
	struct MSData *currPtr;
	
	testMass = (bMass - gCO) / (REAL_4)currentCharge;
	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		if(currPtr->mOverZ > bMass)
		{
			break;
		}
		if(currPtr->mOverZ >= testMass - gParam.fragmentErr &&
			currPtr->mOverZ <= testMass + gParam.fragmentErr)
		{
			test = TRUE;
		}
		currPtr = currPtr->next;
	}
	return(test);
}
/********************************IsThisPossible************************************************
*
*	This function tests to see if the ion is the precursor, precursor - water, below m/z 115,
*	or equal to 120, 136, 147, 159, or 175.  It also checks to make sure that there is sufficient
*	mass to hold multiple charges and is less than the molecular weight of the peptide.
*	currentCharge = current charge.
*	chargeState = charge on the precursor ion.
*/

char IsThisPossible(REAL_4 bMass, INT_4 currentCharge)
{
	REAL_4 minMOverZ, precursor, minWater;
	INT_4	minMassPerCharge = MIN_MASS_PER_CHARGE * gMultiplier + 0.5;
	char maxCharge;
	char test = TRUE;

	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	minMOverZ = (currentCharge - 1) * minMassPerCharge;
		
	if(bMass < minMOverZ)	/*Check that the ion is more than MIN_MASS_PER_CHARGE.*/
	{
		test = FALSE;
	}
	
/* Check to see if the ion is the precursor or the precursor minus water ion.*/
	precursor = (gParam.peptideMW + (maxCharge * gElementMass_x100[HYDROGEN])) / maxCharge;
	
	if(bMass <= (precursor + gParam.fragmentErr) && bMass >= (precursor - gParam.fragmentErr))
	{
		test = FALSE;
	}
	
	minWater = (gParam.peptideMW - gWater + (maxCharge * gElementMass_x100[HYDROGEN])) / maxCharge;
	
	if(bMass <= (minWater + gParam.fragmentErr) && bMass >= (minWater - gParam.fragmentErr))
	{
		test = FALSE;
	}
	
/*	Don't use any ions less than 115 Da.*/
	if(bMass < 115 * gMultiplier)
	{
		test = FALSE;
	}
	
/* Don't use specific ions that are usually immonium ions or tryptic y ions.*/
	if(bMass <= ((gMonoMass_x100[F] - gCO + gElementMass_x100[HYDROGEN]) + gParam.fragmentErr) 
		&& bMass >= ((gMonoMass_x100[F] - gCO + gElementMass_x100[HYDROGEN]) - gParam.fragmentErr))
	{
		test = FALSE;
	}
	
	if(bMass <= (129 * gMultiplier + gParam.fragmentErr) && bMass >= (129 * gMultiplier - gParam.fragmentErr))
	{
		test = FALSE;
	}

	if(bMass <= ((gMonoMass_x100[Y] - gCO + gElementMass_x100[HYDROGEN]) + gParam.fragmentErr) 
		&& bMass >= ((gMonoMass_x100[Y] - gCO + gElementMass_x100[HYDROGEN]) - gParam.fragmentErr))
	{
		test = FALSE;
	}
	
	if(bMass <= ((gMonoMass_x100[W] - gCO + gElementMass_x100[HYDROGEN]) + gParam.fragmentErr) 
		&& bMass >= ((gMonoMass_x100[W] - gCO + gElementMass_x100[HYDROGEN]) - gParam.fragmentErr))
	{
		test = FALSE;
	}
	
	bMass = (bMass * currentCharge) - ((currentCharge - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to singly charged ion.*/
	
	if(bMass > (gParam.peptideMW - gMonoMass_x100[G]))	/*Check that the calculated singly charged b ion is less than 
							the molecular weight of the peptide minus the mass of glycine.*/
	{
		test = FALSE;
	}

	return(test);
}

/*********************************	FindTrypticBIons***************************************************
*
*	This function assumes that the CID ions are all of type b.  The nominal mass values are 
*	determined and the corresponding positions in the array sequenceNodeN are assigned the 
*	additional value of gWeightedIonValues.b.  
*/

void FindTrypticBIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN)
{
	struct MSData *currPtr;
	INT_4 bMassMin, bMassMax, i, j, testForChar;
	REAL_4 bMass, precursor;
	REAL_8 aToMFactor;
	char test, mostLikelyFragCharge, maxCharge;
	
	currPtr = firstMassPtr;
	precursor = gParam.peptideMW / gParam.chargeState;	/*Not exactly the precursor, but close 
														enough.*/

	if(gParam.chargeState == 1)	/*Figure out the most likely charge state of a fragment ion.  This is
							used to determine the gWeightedIonValues.y.*/
	{
		mostLikelyFragCharge = 1;
	}
	else
	{
		mostLikelyFragCharge = gParam.chargeState - 1;
	}
	if(gParam.maxent3)
	{
		mostLikelyFragCharge = 1;
	}
	
	/*max charge state for maxent3 is +1*/
	if(gParam.maxent3)
	{
		maxCharge = 1;
	}
	else
	{
		maxCharge = gParam.chargeState;
	}
	
	while(currPtr != NULL)
	{
		for(i = 1; i <= maxCharge; i++)
		{
			bMass = currPtr->mOverZ;
			test = IsThisPossible(bMass, i);

/*For ions > precursor, check if there is an a ion.*/
			if(test && currPtr->mOverZ > precursor + gParam.fragmentErr * 2)
			{
				test = IsThisStillPossible(bMass, i, firstMassPtr);
			}
			
			if(test)
			{
				bMass = (bMass * i) - ((i - 1) * gElementMass_x100[HYDROGEN]);	/*Convert to +1 ion.*/
/*Alter the values so that they are closer to the expected nominal masses.*/
				if(bMass >= gParam.monoToAv)	/*convert to monoisotopic*/
				{
					aToMFactor = 0;
				}
				else
				{
					if(bMass > (gParam.monoToAv - gAvMonoTransition))
					{
						aToMFactor = (gParam.monoToAv - bMass) / gAvMonoTransition;
					}
					else
					{
						aToMFactor = 1;
					}
				}
				aToMFactor = ((1 - AV_TO_MONO) * aToMFactor) + AV_TO_MONO;
				if(bMass > (gParam.monoToAv - gAvMonoTransition))
				{
					bMass = bMass * aToMFactor;
				}
			
/*The two extremes for the possible b ion nodes are identified first.*/
/*First I'll find the highest mass node.*/
				bMass = bMass + gParam.fragmentErr;
				bMassMax = bMass;		/*Truncate w/o rounding up.*/
	
/*Next I'll find the lowest mass node.*/
				bMass = bMass - gParam.fragmentErr - gParam.fragmentErr + 0.5;
				bMassMin = bMass;	/*Truncate after rounding up.*/

/*Now fill in the middle (if there is anything in the middle.*/
				if(bMassMax >= gGraphLength)
				{
					printf("FindTrypticBIons:  bMassMax >= gGraphLength\n");
					exit(1);
				}
				for(j = bMassMin; j <= bMassMax; j++)
				{
				
					testForChar = sequenceNodeN[j];	/*test to make sure the value fits in a char*/
					
					if(i <= mostLikelyFragCharge)
					{
						if(bMass < precursor && bMass > 147 * gMultiplier)
						{
							testForChar += (INT_4)gWeightedIonValues.b;
						}
						else
						{
							testForChar += (INT_4)(gWeightedIonValues.b * HIGH_MASS_B_MULT);
						}
					}
					else
					{
						if(bMass < precursor && bMass > 147 * gMultiplier)
						{
							testForChar += (INT_4)(gWeightedIonValues.b * HIGH_CHARGE_MULT);
						}
						else
						{
							testForChar += (INT_4)(gWeightedIonValues.b * HIGH_CHARGE_MULT * HIGH_MASS_B_MULT);
						}
					}
					if(testForChar < 127 && testForChar > -127)
					{
						sequenceNodeN[j] = testForChar;
					}
					else
					{
						sequenceNodeN[j] = 63;	/*if 63 is in seqNodeN and C then sum is still <127*/
					}
				}
			}
		}
		currPtr = currPtr->next;
	}
	
/*	If a sequenceNodeN is assigned a non-zero value that is less than the full gWeightedIonValues.b,
*	then that means that it was because an ion was assumed to be of an unlikely charge state.
*	For example, if a doubly charged precursor had an ion that was assumed to be a doubly
*	charged fragment, but the corresponding singly charged ion was absent, then the value
*	for that node would be less than gWeightedIonValues.y.  I remove these from the list here.
*/
	for(i = 0; i < gGraphLength; i++)	
	{
		if(sequenceNodeN[i] < gWeightedIonValues.b)
		{
			sequenceNodeN[i] = 0;
		}
	}


	return;
}

/*********************************TrypticTemplate***********************************************
*
*	This function takes the initialized sequenceNode and modifies it according to empirical 
*	rules developed from years of experience interpretting CID spectra of multiply-charged
*	tryptic peptides.  The rules are as follows:
*		Assume each ion is both a b and a y ion.  Look for an a ion if a b ion is present.
*	Look for b-17 ion if a b ion is present.  Look for an a-17 ion if an a ion is present
*	(and the corresponding b ion is also present).  Look for a y-17 ion if a y ion is present.
*	Look for multiply-charged b and y ions if sufficient mass is available.
*	
*/

void TrypticTemplate(struct MSData *firstMassPtr, SCHAR *sequenceNodeC, 
					SCHAR *sequenceNodeN)
{
	char *ionPresent;
	INT_4 i;
	
	ionPresent = (char *) malloc(gGraphLength * sizeof(char));	/*Will contain C-terminal evidence.*/
	if(ionPresent == NULL)
	{
		printf("TrypticTemplate:  Out of memory");
		exit(1);
	}

	

	
	for(i = 0; i < gGraphLength; i++)	/*Init this array to zero.  These arrays are used to 
										keep track of where the a ions are.*/
	{
		ionPresent[i] = 0;
	}
	
	FindTrypticBIons(firstMassPtr, sequenceNodeN);
				
	FindTrypticB17Ions(firstMassPtr, sequenceNodeN);
						
	FindTrypticAIons(firstMassPtr, sequenceNodeN, ionPresent);
						
	FindTrypticA17Ions(firstMassPtr, sequenceNodeN, ionPresent);
						
	FindTrypticYIons(firstMassPtr, sequenceNodeC);
						
	FindTrypticY17Ions(firstMassPtr, sequenceNodeC);
	
	free(ionPresent);
	
	return;
}


/*******************************SequenceNodeInit************************************************
*
*	This function initializes the array sequenceNode by first assigning the value of zero
*	to each element in the array (from 0 -> 3999).  Next an N-terminal value is assigned to 
*	position number 1, and then the possible C-terminal nodes are found and given an N-terminal
*	value.  There may be more than one C-terminal node, depending on the mass and the error.
*/

void SequenceNodeInit(SCHAR *sequenceNode, SCHAR *sequenceNodeC,
						SCHAR *sequenceNodeN)
{
	INT_4 i, firstNode;
	REAL_4 lastNode;
	REAL_8 aToMFactor;
	INT_4 lastNodeHigh, lastNodeLow;

	for(i = 0; i < gGraphLength; i++)	/*Initialize sequenceNode to zero's.*/
	{
		sequenceNode[i] = 0;
		sequenceNodeC[i] = 0;
		sequenceNodeN[i] = 0;

	}
	
/*	Find the N-terminal node.  This will equal the nominal mass of the N-terminal group R-NH-.*/

	firstNode = gParam.modifiedNTerm;

	sequenceNodeC[firstNode] = N_NODE_VALUE;	/*Give the N-terminal node an arbitrary value.*/
	sequenceNodeN[firstNode] = N_NODE_VALUE;
	sequenceNode[firstNode] = sequenceNodeC[firstNode] + sequenceNodeN[firstNode];
	
/*Figure out what the C-terminal node(s) are.*/

	lastNode = gParam.peptideMW - gParam.modifiedCTerm;
	
/*Alter the values so that they are closer to the expected nominal masses.*/
	if(lastNode >= gParam.monoToAv)	/*convert to monoisotopic*/
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
	if(lastNode > (gParam.monoToAv - gAvMonoTransition))
	{
		lastNode = lastNode * aToMFactor;
	}
	
/*The two extremes for the possible C-terminal nodes are identified first.*/
/*First I'll find the highest mass node.*/
	lastNode = lastNode + gParam.peptideErr;
	lastNodeHigh = lastNode;		/*Truncate w/o rounding up.*/
	if(lastNodeHigh >= gGraphLength)
	{
		printf("SequenceNodeInit:  lastNodeHigh >= gGraphLength\n");
		exit(1);
	}
	sequenceNodeN[lastNodeHigh] = C_NODE_VALUE;
	sequenceNodeC[lastNodeHigh] = C_NODE_VALUE;
	sequenceNode[lastNodeHigh] = sequenceNodeN[lastNodeHigh] + sequenceNodeC[lastNodeHigh];

/*Next I'll find the lowest mass node.*/
	lastNode = lastNode - gParam.peptideErr - gParam.peptideErr;
	lastNodeLow = lastNode;	/*Truncate after rounding up.*/
	sequenceNodeN[lastNodeLow] = C_NODE_VALUE;
	sequenceNodeC[lastNodeLow] = C_NODE_VALUE;
	sequenceNode[lastNodeLow] = sequenceNodeN[lastNodeLow] + sequenceNodeC[lastNodeLow];

/*Now fill in the middle (if there is anything in the middle.*/
	if((lastNodeHigh - lastNodeLow) > 1)
	{
		for(i = lastNodeLow; i <= lastNodeHigh; i++)
		{
			sequenceNodeN[i] = C_NODE_VALUE;
			sequenceNodeC[i] = C_NODE_VALUE;
			sequenceNode[i] = sequenceNodeC[i] + sequenceNodeN[i];
		}
	}

	return;
}


/**********************************MakeSequenceGraph******************************************
*
*	The general idea here is to assume that each observed ion in the linked list "firstMassPtr"
*	is of one of several possibilities - a, b, c, x, y, z+1, y-2, d, v, w, etc..  Making
*	each of these assumptions for each ion, it is possible to mathematically convert the CID
*	data into a "sequence graph".  Here, the sequence graph's index numbers correspond to the
*	nominal mass values of b-type ions.  The value held at each node (or index number) of the 
*	sequence graph is an estimation of the likelihood that a fragmentation is present at that
*	mass.  For example, if there are two ions at 200 and 801, and the peptide molecular weight
*	is 999, then 200 could be a b ion, and 801 would be the corresponding y ion.  The value
*	for the node 200 would reflect the fact that at least two ions suggest that a cleavage 
*	occured at that mass.  
*
*	There are three ways to calculate the node values - a general peptide template, a tryptic
*	peptide template, and a template where there is an arginine present and the precursor
*	has a single charge.  Each template has certain rules that are followed (and are hard-coded)
*	in order to derive the final sequence graph.
*
*	If a sequence tag is known, then this information is used to derive the sequence graph.
*
*	Here is a list of the input variables used by this function:
*	- firstMassPtr has the CID data.
*	- sequenceNode[GRAPH_LENGTH] is the array that eventually contains the final node values.  
*			It is derived from sequenceNodeN and sequenceNodeC.
*	- sequenceNodeN[] is the array that contains the N-terminal fragment info.
*	- sequenceNodeC[] is the array that contains the C-terminal fragment info.
*	- proteolysis indicates if the peptide was cleaved by trypsin, lys-c, or v8.  If so,
*		then the nodes for a C-terminal lysine, arginine, etc. are given higher values.
*		These values will be equal to gWeightedIonValues.y (see below).
*	- fragmentPattern tells whether to use a general, tryptic, or arg+1 template.
*	- modifiedCTerm is needed to convert a putative C-terminal ion into a b ion node.
*	- tagSequence contains a single letter code sequence for the sequence tag.
*	- tagCMass and tagNMass are the terminal unsequenced masses surrounding the sequence tag.
*	- peptideMWPtr is the pointer to the REAL_4 containing the molecular weight of the peptide.
*	- fragmentErr is the CID fragment ion tolerance.
*	- peptideErr is the tolerance for the peptide mass measurement (in Da).
*	- chargeState is the charge state of the precursor ion.
*	- monoToAv is the mass above which average masses are assumed.
* 	- xxxIonVal are the scores used for each ion type.
*	- edmanPresent indicates (TRUE or FALSE) if there is any edman data available.
*	- totalIonVal contains the sum of all of the gWeightedIonValues.
*/

void MakeSequenceGraph(struct MSData *firstMassPtr, SCHAR *sequenceNode, 
                       SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal)
{
	
/*	
*	Initialize the arrays sequenceNode, sequenceNodeC, and sequenceNodeN so that all values
*	are zero, except for the N and C terminal nodes.  This is the only thing that this entire
*	file/function does to the array 'sequenceNode'.  Most of the action occurs to the other two
*	arrays - sequenceNodeN and sequenceNodeC.
*/
	
	SequenceNodeInit(sequenceNode, sequenceNodeC, sequenceNodeN);
						
/*
*	The TrypticTemplate contains the hard-coded rules for interpretting tryptic multiply charged
*	ions (from triple quads).  It uses the various xxxIonVal's to assign values to positions in the arrays called
*	sequenceNodeC and sequenceNodeN, where the indexing of these arrays corresponds to the nominal
*	masses of real and hypothetical (ie, calculated) b type ions.  The hypothetical b type ions
*	derived from those real ions that are assumed to be C-terminal are given node values that are
*	placed in sequenceNodeC.  Likewise, sequenceNodeN contains values for ions that are assumed to 
*	be N-terminal ions.
*/

	if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q')
	{
		TrypticTemplate(firstMassPtr, sequenceNodeC, sequenceNodeN);
	}
	
/*	To be written at a later date.*/

/*	if(fragmentPattern == 'G')
	{
		GeneralTemplate(firstMassPtr, sequenceNode, peptideMW, fragmentErr, chargeState, monoToAv,
						gWeightedIonValues.b, gWeightedIonValues.a, gWeightedIonValues.c, gWeightedIonValues.d, gWeightedIonValues.b_minus17or18, gWeightedIonValues.a_minus17or18, 
						gWeightedIonValues.y, gWeightedIonValues.y_minus2, gWeightedIonValues.y_minus17or18, gWeightedIonValues.x, gWeightedIonValues.z_plus1, gWeightedIonValues.w, 
						gWeightedIonValues.v, gWeightedIonValues.b_minusOH, gWeightedIonValues.b_minusOH_minus17);
	}*/
	
/*
*	The TrypticLCQTemplate contains the hard-coded rules for interpretting tryptic multiply charged
*	ions (from ion traps).  It uses the various xxxIonVal's to assign values to positions in the arrays called
*	sequenceNodeC and sequenceNodeN, where the indexing of these arrays corresponds to the nominal
*	masses of real and hypothetical (ie, calculated) b type ions.  The hypothetical b type ions
*	derived from those real ions that are assumed to be C-terminal are given node values that are
*	placed in sequenceNodeC.  Likewise, sequenceNodeN contains values for ions that are assumed to 
*	be N-terminal ions.
*/

	
	if(gParam.fragmentPattern == 'L')
	{
		TrypticLCQTemplate(firstMassPtr, sequenceNodeC, sequenceNodeN);
	}
/*	
*	RemoveSillyNodes removes all nodes below 260 that cannot be made of any combination of 
*	amino acids.
*/

	RemoveSillyNodes(sequenceNodeC, sequenceNodeN);
	
/*	
*	If a type of proteolysis has been entered by the user (via the char 'proteolysis'), then
*	it is assumed that all possible cleavages characteristic of that proteolysis are present.
*	For example, if the peptide is derived from a tryptic digest, then the node positions 
*	corresponding to the C-terminal node minus 128 (lysine) and 156 (arginine) are calculated.
*	if the values in sequenceNodeN[C-term minus residue] and sequenceNodeC[C-term minus residue]
*	are both zero (ie, no real ions have indicated a cleavage at that site), then these arrays
*	are assigned a value of one at that node value.  If the data is of poor quality, this permits
*	the possibility of continuing a sequence to the C-terminus even if the cleavage info is
*	absent.  By assigning a value of one to those nodes, I do not weight these sequences very 
*	heavily.
*/
	
	if(gParam.proteolysis != 'N')
	{
		AddCTermResidue(sequenceNodeC, sequenceNodeN);
	}
	
/*	
*	If 'edmanPresent' is TRUE and the N-terminus is not modified, then call the function 
*	AddEdmanData.  AddEdmanData uses the globals gMaxCycleNum and 
*	gEdmanData[MAX_PEPTIDE_LENGTH][gAminoAcidNumber] to alter the arrays 'sequenceNodeC' 
*	and 'sequenceNodeN'.  
*/

	if(gParam.edmanPresent)
	{
		AddEdmanData(sequenceNodeC, sequenceNodeN, totalIonVal);
	}

/*	
*	If a sequence tag has been provided, then the region containing the sequence is excised
*	and the graph closed over the wound.  The peptideMW loses the mass of the sequence tag.
*/
	
	if(gParam.tagSequence[0] != '*')
	{
		AddTag(sequenceNodeC, sequenceNodeN);
	}

	return;
}
