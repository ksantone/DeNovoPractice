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
	
	LutefiskSummedNode is a file containing the function SummedNodeScore plus its associated
	functions.  It was written to be used as part of a program called "LutefiskXP", which
	is used to aid in the interpretation of CID data of peptides.  The general aim of this 
	(and the function SummedNodeScore) is to connect the nodes starting from the C-terminal
	node(s), where the differences between the nodes corresponds to the nominal mass of an
	amino acid residue.  If a node that is connected to the C-terminus is found that cannot
	be extended any further towards the N-terminus, then it is saved as a "one-edged node".  
	These special nodes are used later for making two amino acid jumps (for those situations
	where there is no fragmentation at particular peptide bonds).
*/

#include <stdio.h>
#include <stdlib.h>
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"

/*Define some globals that were defined in LutefiskMain.c.*/
extern INT_4 gNomMass[AMINO_ACID_NUMBER];

/********************************AddExtraNodes***********************************************
*
*	This function compares the arrays "sequenceNode" and "evidence" to see if there are
*	nodes that were not connected to the C-terminus.  It ignores nodes that are adjacent
*	to nodes that had been connected.  It adds the values in sequenceNodeN and sequenceNodeC
*	and puts them into the array sequenceNode.
*/
void AddExtraNodes(SCHAR *sequenceNode, SCHAR *sequenceNodeN, 
		   SCHAR *sequenceNodeC, char *evidence)
{
	INT_4 i, j, k, newNodeValue;
	char test;
	
	i = gGraphLength - 1;	/*Start out at the maximum node position.*/
	while(i >= 0)	/*Don't let the index value go negative (there are no negative b ions).*/
	{
		if(evidence[i] != 0)	/*If I hit a node that has any cleavage evidence, then do this.*/
		{
			test = TRUE;	/*Test is used to determine if consecutive evidence nodes have been
							represented already in the sequenceNode array.  In other words, if
							there is a series of three consecutive nodes w/ non-zero evidence,
							it may be that the middle one actually can be connected to the
							C-terminus (as above), and therefore that node position or index 
							value for the array sequenceNode is non-zero.  In such situations,
							I don't bother adding the adjacent nodes to the array sequenceNode.*/
			j = i;	/*I need a new index to step through the consecutive non-zero nodes.*/
			while(evidence[j] != 0)
			{
				if(sequenceNode[j] != 0)
				{
					test = FALSE;	/*I've found a node that has already been included in sequenceNode.*/
				}
				j--;	/*Keep stepping down so that I can find the end of this series of nodes.*/
			}
			if(test || gParam.fragmentErr <= 0.5 * gMultiplier)	/*Alter the relevant values of sequenceNode.
													If the fragment tolerance is less than 0.5,
													then adjacent nodes are not due to slop in
													mass assignments, and should therefore be
													taken seriously.*/
			{
				for(k = j + 1; k <= i; k++)	/*evidence[j] is zero, so just go up to it.*/
				{
					newNodeValue = sequenceNodeC[k] + sequenceNodeN[k];
					
					if(newNodeValue > 127 || newNodeValue < -127)
					{
						newNodeValue = 127;
					}
					
					if(newNodeValue > sequenceNode[k])
					{
						sequenceNode[k] = newNodeValue;
					}
				}
			}
					
			i = j + 1;	/*Must be j+1, cuz its incremented down before starting the while loop.*/
		}
		i--;
	}
	
	
	return;
}

/********************************SortOneEdgeNodes*******************************************
*
*	This function sorts the list of one edge nodes and removes the redundancies and sorts by
*	increasing mass.
*/

void SortOneEdgeNodes(INT_4 *oneEdgeNodes, INT_4 *oneEdgeNodesIndex)
{
	INT_4 i;
	char test = TRUE;
	INT_4 *tempNode, tempIndex, smallestTemp, smallestTempValue;
	
	tempNode = (int *) malloc(gGraphLength * sizeof(INT_4));	/*Will contain C-terminal evidence.*/
	if(tempNode == NULL)
	{
		printf("SortOneEdgeNodes:  Out of memory");
		exit(1);
	}
	
	if(*oneEdgeNodesIndex >= gGraphLength)
	{
		printf("SortOneEdgeNodes:  *oneEdgeNodesIndex >= gGraphLength\n");
		exit(1);
	}

	
	tempIndex = 0;
	
	while(test)
	{
		test = FALSE;
		
		i = 0;
		while(i < *oneEdgeNodesIndex)	/*Find a positive number.*/
		{
			if(oneEdgeNodes[i] > 0)
			{
				test = TRUE;
				smallestTemp = i;
				smallestTempValue = oneEdgeNodes[i];
				break;
			}
			i++;
		}
		
		if(test)
		{		
			for(i = 0; i < *oneEdgeNodesIndex; i++)	/*Find the smallest value that is not zero.*/
			{
				if(oneEdgeNodes[i] > 0)
				{
					if(oneEdgeNodes[i] < smallestTempValue)
					{
						smallestTemp = i;
						smallestTempValue = oneEdgeNodes[i];
					}
				}
			}
			
			for(i = 0; i < *oneEdgeNodesIndex; i++)	/*Find the ones that are identical.*/
			{
				if(i != smallestTemp)
				{
					if(smallestTempValue == oneEdgeNodes[i])
					{
						oneEdgeNodes[i] = 0;
					}
				}
			}
			
			tempNode[tempIndex] = smallestTempValue;	/*Set the temporary value.*/
			tempIndex++;							/*Increment the number of values.*/
			if(tempIndex >= gGraphLength)
			{
				printf("gGraphLength is too small.");
				exit(1);
			}
			oneEdgeNodes[smallestTemp] = 0;			/*Set this to zero so that its not used
													again.*/
		}
	}
	
	*oneEdgeNodesIndex = tempIndex;	/*Transfer to the oneEdgeNodes array.*/
	if(*oneEdgeNodesIndex >= gGraphLength)
	{
		printf("SortOneEdgeNodes:  *oneEdgeNodesIndex >= gGraphLength\n");
		exit(1);
	}
	for(i = 0; i < tempIndex; i++)
	{
		oneEdgeNodes[i] = tempNode[i];
	}

	free(tempNode);

	return;
}

/********************************FindCurrentNode*********************************************
*
*	This function takes the currentNode from which amino acid residue masses have already been
*	subtracted, and finds the next node lower in mass that can be connected to the C-terminus.
*	That new node serves as the next currentNode.
*/

INT_4 FindCurrentNode(SCHAR *sequenceNode, INT_4 currentNode)
{
	INT_4 i;
	
	i = currentNode;
	
	while(i > 0)
	{
		i--;
		
		if(i < currentNode - 190 * gMultiplier)	/* If i goes below the mass of currentNode - Trp, then
									terminate the search.  I used 190 as the mass of Trp 
									plus an inordinately large error.*/
		{
			currentNode = 0;
			return(currentNode);
		}
		
		if(sequenceNode[i] > 0)
		{
			currentNode = i;
			return(currentNode);
		}
	}
	
	currentNode = i;

	return(currentNode);
}

/********************************AssignProNodeValue*********************************************
*
*	This function figures out how to score the node.  If this function has been called,
*	then its already known that this is one of the connectable nodes, but the score to
*	be placed in the array sequenceNode now needs to be determined.  This is based on the
*	type of evidence for the nextNode and the currentNode.  If the nextNode evidence is 
*	that both N- and C-terminal ions were present at that point, then the values from
*	sequenceNodeC,sequenceNodeN, and sequenceNode[currentNOde] are summed.  If the
*	evidence for the nextNode does not match w/ the currentNode, then the values for 
*	sequenceNodeC and sequenceNodeN only are summed.
*/
void AssignProNodeValue(INT_4 nextNode, INT_4 currentNode, char *evidence, 
		        SCHAR *sequenceNode, SCHAR *sequenceNodeC,
		  	SCHAR *sequenceNodeN, INT_4 totalIonVal)
{
	INT_4 nodeScore = 0;
		
/*	The 2 aa extension w/ proline is only used if its a tryptic peptide and therefore
	only y ions would be expected if there are no ions delineating pro.*/
	if((evidence[nextNode] == 'B' || evidence[nextNode] == 'C') &&
		(evidence[currentNode] == 'B' || evidence[currentNode] == 'C'))
	{
		nodeScore = sequenceNodeC[nextNode] + sequenceNodeN[nextNode]
					+ (totalIonVal * TOTALIONVAL_MULTIPLIER * 0.5);
					
	}

	if(sequenceNode[nextNode] < 0)	/*If the nextNode value is negative, make it positive.*/
	{
		sequenceNode[nextNode] = -1 * sequenceNode[nextNode];
	}
	
	if(nodeScore > 127 || nodeScore < -127)
	{
		nodeScore = 127;
	}
	
	if(nodeScore > sequenceNode[nextNode])	/*If the new nodeScore is greater, then assign its
											value to sequencNode[nextNode], otherwise leave
											the original value as a positive number.*/
	{
		sequenceNode[nextNode] = nodeScore;
	}
		
	return;
}


/********************************AssignNodeValue*********************************************
*
*	This function figures out how to score the node.  If this function has been called,
*	then its already known that this is one of the connectable nodes, but the score to
*	be placed in the array sequenceNode now needs to be determined.  This is based on the
*	type of evidence for the nextNode and the currentNode.  If the nextNode evidence is 
*	that both N- and C-terminal ions were present at that point, then the values from
*	sequenceNodeC,sequenceNodeN, and sequenceNode[currentNOde] are summed.  If the
*	evidence for the nextNode does not match w/ the currentNode, then the values for 
*	sequenceNodeC and sequenceNodeN only are summed.
*/
void AssignNodeValue(INT_4 nextNode, INT_4 currentNode, char *evidence, 
					SCHAR *sequenceNode, SCHAR *sequenceNodeC,
					SCHAR *sequenceNodeN, INT_4 totalIonVal)
{
	INT_4 nodeScore = 0;
	INT_4 i = 0;
	REAL_4 scoreAdjuster;
	
	/*if(nextNode >= 10428 && nextNode <= 10434)
	{
		i++;	
		i++;
	}	for debugging*/
	
	scoreAdjuster = currentNode - nextNode; /*Nominal amino acid mass of this connection.*/
	scoreAdjuster = scoreAdjuster / gAvResidueMass; /*Ratio between this connection and
														the average connection.*/
	scoreAdjuster = (scoreAdjuster + 99) / 100;	/*Reduce scoreAdjuster effect.*/
	
	if(evidence[nextNode] == 'B' || evidence[currentNode] == 'B')
	{
		nodeScore = sequenceNodeC[nextNode] + sequenceNodeN[nextNode]
					+ (totalIonVal * TOTALIONVAL_MULTIPLIER);
					/*((nodeScore + totalIonVal) + (nodeScore * totalIonVal)) / 2;*/
					/*(totalIonVal * TOTALIONVAL_MULTIPLIER);*/
					/*+ sequenceNode[currentNode] * TOTALIONVAL_MULTIPLIER;*/
	}
	else
	{
		if(evidence[nextNode] != evidence[currentNode])
		{
			nodeScore = sequenceNodeC[nextNode] + sequenceNodeN[nextNode];
		}
		else
		{
			nodeScore = sequenceNodeC[nextNode] + sequenceNodeN[nextNode]
						+ (totalIonVal * TOTALIONVAL_MULTIPLIER);
						/*((nodeScore + totalIonVal) + (nodeScore * totalIonVal)) / 2;*/
						/*(totalIonVal * TOTALIONVAL_MULTIPLIER);*/
						/*+ sequenceNode[currentNode] * TOTALIONVAL_MULTIPLIER;*/
		}
	}
	
	nodeScore = (nodeScore * scoreAdjuster) + 0.5;	/*Adjust score for size of extension.*/
	
	if(sequenceNode[nextNode] < 0)	/*If the nextNode value is negative, make it positive.*/
	{
		sequenceNode[nextNode] = -1 * sequenceNode[nextNode];
	}
	
	if(nodeScore > 127 || nodeScore < -127)
	{
		nodeScore = 127;
	}
	
	if(nodeScore > sequenceNode[nextNode])	/*If the new nodeScore is greater, then assign its
											value to sequencNode[nextNode], otherwise leave
											the original value as a positive number.*/
	{
		sequenceNode[nextNode] = nodeScore;
	}
		
	return;
}

/***************************************InitSummedNodeArrays********************************
*
*		This function initializes some of the arrrays used by the function SummedNodeScore.
*/
void InitSummedNodeArrays(SCHAR *sequenceNode, SCHAR *sequenceNodeC, 
							SCHAR *sequenceNodeN, INT_4 *oneEdgeNodes, 
							char *evidence)
{
	INT_4 i;
	
	for(i = 0; i < gGraphLength; i++)	/*Initialize these arrays.*/
	{
		sequenceNode[i] = 0;
		evidence[i] = 0;
	}
	
	for(i = 0; i < gGraphLength; i++)
	{
		oneEdgeNodes[i] = 0;
	}
	
	for(i = 0; i < gGraphLength; i++)	/*Set up the evidence array.*/
	{
	
		if(sequenceNodeC[i] != 0 && sequenceNodeN[i] != 0)
		{
			evidence[i] = 'B';
		}
	
		if(sequenceNodeC[i] != 0 && sequenceNodeN[i] == 0)
		{
			evidence[i] = 'C';
		}
		
		if(sequenceNodeC[i] == 0 && sequenceNodeN[i] != 0)
		{
			evidence[i] = 'N';
		}
	}
	
	for(i = 0; i < gGraphLength; i++)	/*-1 implies a sequencetag*/
	{
		if(sequenceNodeC[i] == -1 && sequenceNodeN[i] == -1)
		{
			evidence[i] = 'B';
		}
	}
	return;
}

	
/******************************SummedNodeScore**********************************************
*
*		This function is used to connect the nodes that were identified in the function
*	MakeSequenceGraph.  The input data is stored in the INT_4 arrays sequenceNodeC
*	and sequenceNodeN.  The array sequenceNodeC contains all of the evidence for C-terminal
*	ions, and the array sequenceNodeN contains the evidence for N-terminal ions.  Both
*	arrays were formed by assuming that each observed ion could be one of several possiblities -
*	b, a, y, etc. - and then mathematically converting to the corresponding b ion value. 
*		Here the program starts to connect the nodes, beginning at the C-terminus.  It does
*	not remember exactly how it connects the nodes, rather the point is to assign high score
*	values to the array sequenceNode for those nodes where its possible to connect with the
*	C-terminus.  Later, the program will start making subsequences starting at the N-terminus,
*	and it will use the scores held in sequenceNode to guide the way towards the C-terminus.
*		In addition, this function keeps track of those nodes that can be connected to the
*	C-terminus, but do not connect to any other nodes N-terminal to it.  These so-called one-
*	edged nodes will be used for making two amino acid jumps between nodes.
*/

void SummedNodeScore(SCHAR *sequenceNode, SCHAR *sequenceNodeC, 
					SCHAR *sequenceNodeN, INT_4 *oneEdgeNodes,
					INT_4 *oneEdgeNodesIndex, INT_4 totalIonVal)
{
	INT_4 i, k;
	INT_4 j, currentNode, nextNode, lowSuperNode, highSuperNode;
	INT_4 highArgNode, highLysNode, lowArgNode, lowLysNode, lowMassCutoff, y2;
	char *evidence;
	char test = TRUE;
	char anotherTest, doIt;	

	
	evidence = (char *) malloc(gGraphLength * sizeof(char ));	/*Will contain C-terminal evidence.*/
	if(evidence == NULL)
	{
		printf("SummedNodeScore:  Out of memory");
		exit(1);
	}
	
	*oneEdgeNodesIndex = 0;
	
/*	Find the low mass cutoff to be used for LCQ data.*/
	lowMassCutoff = (gParam.peptideMW + (gParam.chargeState * gElementMass_x100[HYDROGEN])) / 
					gParam.chargeState;
	lowMassCutoff = lowMassCutoff * 0.333;	/*Use the rule of 1/3 for the low mass end cutoff.*/
	
/*	Find the lowSuperNode and highSuperNode*/
	highSuperNode = gGraphLength -1;		/*default values when no specific sequencetag is used*/
	lowSuperNode = 0;
	
	if(gParam.tagSequence[0] != '*')
	{
		for(i = gGraphLength - 1; i >= 0; i--)
		{
			if(sequenceNodeN[i] == -1 && sequenceNodeC[i] == -1)
			{
				highSuperNode = i;
				break;
			}
		}
		
		for( i = 0; i < gGraphLength; i++)
		{
			if(sequenceNodeN[i] == -1 && sequenceNodeC[i] == -1)
			{
				lowSuperNode = i;
				break;
			}
		}
	}
	
/*	
*	Find the nodes corresponding to the C-terminal node minus R or K.  This is used for
*	a section below that compensates for the absence of y2 ions in higher mass peptides
*	obtained from ion traps w/ a low mass cutoff of 1/3 of the precursor ion.
*/
	i = gGraphLength - 1;
	while(i > 0 && sequenceNodeC[i] == 0)
	{
		i--;
	}
	highArgNode = i - gGapList[R];
	highLysNode = i - gGapList[K];
	while(i > 0 && sequenceNodeC[i] != 0)
	{
		i--;
	}
	i = i + 1;	/*Go back up to where evidence is not zero.*/
	lowArgNode = i - gGapList[R];
	lowLysNode = i - gGapList[K];
	
	
/*	Initialize sequenceNode, oneEdgeNodes, and evidence to zero.  Then assign values of
*	either B, N, or C to those nodes that have any evidence derived from both the N- and
*	C-terminii, the N-terminus only, or the C-terminus only.
*/
	InitSummedNodeArrays(sequenceNode, sequenceNodeC, sequenceNodeN, oneEdgeNodes, evidence);

	
/*	Make sure that the index is greater than zero.  The first time evidence is not zero,
*	the value of test becomes FALSE, which means that thereafter the while loop will continue
*	only as INT_4 as evidence does not equal zero, which is only possible for those nodes
*	that are due to multiple C-terminal nodes resulting from a large error in the peptide
*	mass measurement.
*/
	i = gGraphLength - 1;
	while(i > 0 && (test || evidence[i] != 0))
	{
		i--;
		if(evidence[i] != 0)
		{
			test = FALSE;
			currentNode = i;

/*	Assign the value to the C-terminal node of sequenceNode.*/
			if(sequenceNodeC[currentNode] + sequenceNodeN[currentNode] < 127)
			{
				sequenceNode[currentNode] = sequenceNodeC[currentNode] + sequenceNodeN[currentNode];
			}
			else
			{
				sequenceNode[currentNode] = 127;
			}

/*	Now that I've found a C-terminal Node, lets start connecting the dots.*/

/*
*	A connection is made
*	when the difference between the currentNode and a lower value node equals the nominal mass
*	of an amino acid residue.  That lower value node must have evidence of being a real cleavage
*	in that the value of evidence[node] is not zero (ie, is B, C, or N).  Once a connection
*	is made, then a value for that node is assigned to the array sequenceNode[node].  If a 
*	previous connection has made an assignment to the array sequenceNode at the node under 
*	consideration, then the node with the greatest absolute value is kept (and is made positive,
*	since it is now part of the current set of node connections}.  Once all
*	of the connections have been made and sequenceNode assignments have been completed, then
*	all values are given a negative value in order to distinguish old node connections with
*	any new ones using differnt C-terminii.  For a given currentNode, the program checks the
*	twenty possible nextNode values.  After that it finds the next currentNode, which will be
*	the next node down that has a positive value for sequenceNode[node].
*/

			while(currentNode != 0)
			{
				anotherTest = TRUE;
				for(j = 0; j < gAminoAcidNumber; j++)
				{
					if(gGapList[j] != 0)
					{
						nextNode = currentNode - gGapList[j];
						doIt = TRUE;
						if(currentNode > highSuperNode && nextNode < lowSuperNode)
						{
							doIt = FALSE;	/*you skipped a superNode*/
						}
						if(nextNode >= 0 && evidence[nextNode] != 0 && doIt)
						{
							anotherTest = FALSE;
							AssignNodeValue(nextNode, currentNode, evidence, sequenceNode, 
									sequenceNodeC, sequenceNodeN, totalIonVal);
						}
					}
				}
				
				if(gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q' 
					|| gParam.fragmentPattern == 'L')
				{
					for(j = 0; j < gAminoAcidNumber; j++)	/*check for 2aa's w/ proline*/
					{
						if(gGapList[j] != 0)
						{
							nextNode = currentNode - gGapList[j] - gGapList[P];
							doIt = TRUE;
							if(currentNode > highSuperNode && nextNode < lowSuperNode)
							{
								doIt = FALSE;	/*you skipped a superNode*/
							}
							if(nextNode >= 0 && evidence[nextNode] != 0 && doIt)
							{
								anotherTest = FALSE;
								AssignProNodeValue(nextNode, currentNode, evidence, sequenceNode, 
												sequenceNodeC, sequenceNodeN, totalIonVal);
							}
						}
					}
				}
				
				/*For LCQ data the y2 ions may be missing for ions greater than 1200 
				(1200 is chosen because GK y2 would be below the ion trap low mass cutoff),
				so 2 amino acid extensions are allowed for tryptic peptides below R or K.*/
				if(gParam.fragmentPattern == 'L' && gParam.peptideMW > 1200 * gMultiplier
					&& gParam.proteolysis == 'T' && gParam.chargeState <= 2)
				{
					if((currentNode <= highArgNode && currentNode >= lowArgNode) ||
						(currentNode <= highLysNode && currentNode >= lowLysNode))
					{
						for(j = 0; j < gAminoAcidNumber; j++)	/*check for 2aa's*/
						{
							if(gGapList[j] != 0)
							{
								for(k = j; k < gAminoAcidNumber; k++)
								{
									if(gGapList[k] != 0)
									{
										nextNode = currentNode - gGapList[j] - gGapList[k];
										doIt = TRUE;
										if(currentNode > highSuperNode && nextNode < lowSuperNode)
										{
											doIt = FALSE;	/*you skipped a superNode*/
										}
										y2 = 147 * gMultiplier;	/*y1 ion for c-term Lys*/
										if(gGapList[j] < gGapList[k])
										{
											y2 += gGapList[j];	/*add the lowest mass aa*/
										}
										else
										{
											y2 += gGapList[k];
										}
										if(y2 > lowMassCutoff)
										{
											doIt = FALSE;	/*this y2 should be within range*/
										}
										if(nextNode >= 0 && evidence[nextNode] != 0 && doIt)
										{
											anotherTest = FALSE;
											AssignNodeValue(nextNode, currentNode, evidence, 
																sequenceNode, sequenceNodeC, 
																sequenceNodeN, totalIonVal);
										}
									}
								}
							}
						}
					}
				}


				if(anotherTest)
				{
					if(*oneEdgeNodesIndex < gGraphLength - 1)
					{
						oneEdgeNodes[*oneEdgeNodesIndex] = currentNode;
						*oneEdgeNodesIndex += 1;
						if(*oneEdgeNodesIndex >= gGraphLength)
						{
							printf("SummedNodeScore:  *oneEdgeNodesIndex >= gGraphLength\n");
							exit(1);
						}
					}
					else
					{
						printf("gGraphLength is too small.");
						exit(1);
					}
				}
/*	
*	If I reach the N-terminus, or no further connections can be made, then FindCurrentNode
*	returns a value of zero, which terminates the while loop.
*/
				currentNode = FindCurrentNode(sequenceNode, currentNode);
				if(currentNode > gGraphLength || currentNode < 0)
				{
					printf("SummedNodeScore:  currentNode > gGraphLength || currentNode < 0\n");
					exit(1);
				}
			}
			
			for(j = 0; j < gGraphLength; j++)	/*Make all of the positive values negative.*/
			{
				if(sequenceNode[j] > 0)
				{
					sequenceNode[j] = -1 * (INT_2)sequenceNode[j];
				}
			}
		}
	}
	
	for(i = 0; i < gGraphLength; i++)	/*Make everything positive.*/
	{
		if(sequenceNode[i] < 0)
		{
			sequenceNode[i] = -1 * sequenceNode[i];
		}
	}
	
	SortOneEdgeNodes(oneEdgeNodes, oneEdgeNodesIndex);
	
	AddExtraNodes(sequenceNode, sequenceNodeN, sequenceNodeC, evidence);
	
/*	Add -1 to the superNode positions of sequenceNode*/
	for(i = 0; i < gGraphLength; i++)
	{
		if(sequenceNodeN[i] == -1 && sequenceNodeC[i] == -1)
		{
			sequenceNode[i] = -1;
		}
	}
	
	if(gParam.fMonitor && gCorrectMass)
	{	
		printf("Graph is finished. \n");
	}
	
	free(evidence);	
	
	return;
}
	
	
	
	
	
	
	
	
