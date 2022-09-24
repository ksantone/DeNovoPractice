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

#ifndef __LUTEFISK__
#define __LUTEFISK__

#include "LutefiskDefinitions.h"


/*	Prototypes for LutefiskMain.	*/
void 			Run(void);
void 			ReadParamsFile(void);
INT_4 			ReadDetailsFile(void);
void 			SetupGapList();
void 			ReadEdmanFile();
void 			FreeSequenceScore(struct SequenceScore *currPtr);
void 			FreeMassList(struct MSData *currPtr);
void 			FreeSequence(struct Sequence *currPtr);
void 			CreateGlobalIntegerMassArrays(struct MSData *firstMassPtr);
void 			ChangeOutputName(void);
void 			FindTheMultiplier(void);
static void 	BuildPgmState(int argc, char **argv);
BOOLEAN 		SystemCheck(void);
void	 		ReadResidueFile(void);
void 			SetupSequenceTag();
void 			PrintPartingGiftToFile(void);
void 			PrintHeaderToFile(void);
void			AdjustPeptideMW(struct MSData *firstMassPtr);

/*Prototypes for Haggis*/
struct Sequence *Haggis(struct Sequence *firstSequencePtr , struct MSData *firstMassPtr);
BOOLEAN			NodeStep(INT_4 *nodeNum, INT_4 *nodeMass);
void			StoreSeq(INT_4 nodeNum, INT_4 *nodeMass);
INT_4			ResidueMass(INT_4 inputMass);
INT_4			*LoadMassArrays(INT_4 *mass, struct MSData *firstMassPtr, INT_4 charge);
void			SetupBackwardAndForwardNodes(INT_4 *mass);
void			FindNodeSequences(INT_4 *mass);
void			GetSequenceOfResidues(INT_4 *mass);
struct Sequence	*LinkHaggisSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr);
struct Sequence *LoadHaggisSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
						INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection);
void 			FleshOutSequenceEnds(struct MSData *firstMassPtr);
void			GetCTermMasses(INT_4 *cTermMassNum, INT_4 *cTermMasses);
void			GetNTermMasses(INT_4 *nTermMassNum, INT_4 *nTermMasses);
void			GetBestNtermSeq(INT_4 *sequenceToAdd, INT_4 mass, struct MSData *firstMassPtr);
INT_4			RatchetHaggis(INT_4 *sequence, INT_4 residueNum, INT_4 position, INT_4 ratchetMass, INT_4 correctMass);
void			AppendSequences(void);
void			MakeAAArray(void);
REAL_4			NSequenceScore(INT_4 mass, INT_4 *sequence, INT_4 residueNum, struct MSData *firstMassPtr);
void			GetBestCtermSeq(INT_4 *sequenceToAdd, INT_4 mass, char cTerm, struct MSData *firstMassPtr);
char			CheckCterm(struct MSData *firstMassPtr);
REAL_4			CSequenceScore(INT_4 mass, INT_4 *sequence, INT_4 residueNum, struct MSData *firstMassPtr);
void			ModifyHaggisSequences(struct MSData *firstMassPtr);
BOOLEAN			FindBIon(INT_4 sequenceIndex,INT_4 residueIndex, struct MSData *firstMassPtr);
BOOLEAN			FindYIon(INT_4 sequenceIndex,INT_4 residueIndex, struct MSData *firstMassPtr);


/*	Prototypes for LutefiskGetCID.*/
tMSDataList		*ReadCIDFile(char *inFilename);
void 			AddTheIonOffset(tMSDataList *inMSDataList);
void 			CalibrationCorrection(tMSDataList *inPeakList);
void 			CentroidOrProfile(tMSDataList *inMSDataList);
void 			CheckConnections(tMSDataList *inPeakList);
void 			CheckSignalToNoise(tMSDataList *inPeakList, tMSDataList *inMSDataList);
void 			CheckTheIntensity(tMSDataList *inMSDataList);
void 			DefectCorrection(tMSDataList *inPeakList);
void 			FindTheGoldenBoys(tMSDataList *inMSDataList);
INT_4 			FindThreshold(tMSDataList *inMSDataList);
REAL_4 			GetPeakWidth(tMSDataList *inMSDataList);
void 			GuessAtTheFragmentPattern();
tMSDataList		*IonCondenser(tMSDataList *inMSDataList);
INT_4 			IntensityDescendSortFunc(const void *n1, const void *n2);
INT_4 			MassAscendSortFunc(const void *n1, const void *n2);
void 			NormalizeIntensity(tMSDataList *inMSDataList);
INT_4 			PurgeTheWindow(tMSDataList *inMSDataList, tMSData *windowStartPtr, 
						 INT_4 ionsInWindow, INT_4 endingMass/*, INT_4 maxIonsPerWindow*/);
void 			RemoveIsotopes(tMSDataList *inMSDataList);
void 			RemovePrecursors(tMSDataList *inMSDataList);
void 			SmoothCID(tMSDataList *inMSDataList);
void 			WeedTheIons(tMSDataList *inMSDataList, INT_4 finalIonCount, BOOLEAN spareGoldenBoys);
void 			WindowFilter(tMSDataList *inMSDataList);
void 			FreeAllMSData(struct MSData *currPtr);
void 			free_all(struct MSData *s);
struct MSData 	*LoadMSDataStruct(REAL_4 massValue, INT_4 ionIntensity);
struct MSData 	*AddToCIDList(struct MSData *firstPtr, struct MSData *currPtr);
void 			ModifyList(struct MSData *firstPtr, struct MSData *currPtr);
INT_4 			FindMedian(struct MSData *firstPtr);
char 			*my_fgets(char *s, INT_4 n, FILE *fp);						 
struct MSData 	*GetCidData(void);
struct MSData 	*AddToListNoNull(struct MSData *firstPtr, struct MSData *currPtr);
void 			SortByMass(struct MSData *firstAvMassPtr);
INT_4 			countIons(struct MSData *firstAvMassPtr);
struct MSData 	*ZeroTheIons(struct MSData *firstAvMassPtr);
INT_4 			countIonsAgain(struct MSData *firstAvMassPtr);
void 			EliminateBadHighMassIons(tMSDataList *inPeakList);
void			LowMassIonRemoval(tMSDataList *peakList);
REAL_4 			GeneralEval(tMSDataList  *peakList);
REAL_4 			LowMassIonCheck(tMSDataList  *peakList);
REAL_4 			HighMassIonCheck(tMSDataList  *peakList);
void 			FindBYGoldenBoys(tMSDataList *inMSDataList);

#ifdef DEBUG
void 			DumpMassList(struct MSData *firstPtr);
void 			DumpMSData(tMSDataList *inList);
#endif
					
/*Prototypes for LutefiskMakeGraph.*/
void 			MakeSequenceGraph(struct MSData *firstMassPtr, SCHAR *sequenceNode, 
						SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal);
void 			SequenceNodeInit(SCHAR *sequenceNode, SCHAR *sequenceNodeC, 
									SCHAR *sequenceNodeN);
void 			TrypticTemplate(struct MSData *firstMassPtr, SCHAR *sequenceNodeC, 
						SCHAR *sequenceNodeN);
void 			FindTrypticBIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN);
char 			IsThisPossible(REAL_4 bMass, INT_4 currentCharge);
void 			FindTrypticB17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN);
void 			FindTrypticAIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent);
void 			FindTrypticA17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent);
void 			FindTrypticYIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeC);
void 			FindTrypticY17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeC);
void 			AddTag(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN);
void 			AddCTermResidue(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN);
char 			RatchetIt(INT_4 *aaNum, char cycle, char *sequence, INT_4 seqLength);
void 			AddEdmanData(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal);
char 			IsThisStillPossible(REAL_4 bMass, INT_4 currentCharge, struct MSData *firstMassPtr);
void 			TrypticLCQTemplate(struct MSData *firstMassPtr, SCHAR *sequenceNodeC, 
						SCHAR *sequenceNodeN);
void 			FindTrypticLCQBIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN);
void 			FindTrypticLCQB17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN);
void 			FindTrypticLCQAIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent);
void 			FindTrypticLCQA17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeN, char *ionPresent);
void 			FindTrypticLCQYIons(struct MSData *firstMassPtr, SCHAR *sequenceNodeC);
void 			FindTrypticLCQY17Ions(struct MSData *firstMassPtr, SCHAR *sequenceNodeC);
void 			RemoveSillyNodes(SCHAR *sequenceNodeC, SCHAR *sequenceNodeN);

						
/*Prototypes for LutefiskSummedNode.*/
void 			SummedNodeScore(SCHAR *sequenceNode, SCHAR *sequenceNodeC, 
						SCHAR *sequenceNodeN, INT_4 *oneEdgeNodes, INT_4 *oneEdgeNodesIndex, INT_4 totalIonVal);
void 			InitSummedNodeArrays(SCHAR *sequenceNode, SCHAR *sequenceNodeC, 
						SCHAR *sequenceNodeN, INT_4 *oneEdgeNodes, char *evidence);
void 			AssignNodeValue(INT_4 nextNode, INT_4 currentNode, char *evidence, 
						SCHAR *sequenceNode, SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal);
INT_4 			FindCurrentNode(SCHAR *sequenceNode, INT_4 currentNode);
void 			SortOneEdgeNodes(INT_4 *oneEdgeNodes, INT_4 *oneEdgeNodesIndex);
void 			AddExtraNodes(SCHAR *sequenceNode, SCHAR *sequenceNodeN, SCHAR *sequenceNodeC, char *evidence);
void 			AssignProNodeValue(INT_4 nextNode, INT_4 currentNode, char *evidence, SCHAR *sequenceNode, 
						SCHAR *sequenceNodeC, SCHAR *sequenceNodeN, INT_4 totalIonVal);

/*Prototypes for SubsequenceMaker.*/
struct Sequence *SubsequenceMaker(INT_4 *oneEdgeNodes, INT_4 oneEdgeNodesIndex, 
						SCHAR *sequenceNode);
struct Sequence *NterminalSubsequences(SCHAR *sequenceNode, INT_4 maxLastNode, INT_4 lowSuperNode,
						INT_4 highSuperNode);
struct Sequence *LoadSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
						INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection);
struct Sequence *LoadFinalSequenceStruct(INT_4 *peptide, INT_4 peptideLength, 
						INT_4 score, INT_4 nodeValue, INT_4 gapNum, INT_2 nodeCorrection);
struct Sequence *LinkSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr);
struct Sequence *StoreSubsequences(struct Sequence *newSubsequencePtr, struct extension *extensionList,
						struct Sequence *currentSubsequence,INT_4 *lastNode, INT_4 lastNodeNum,
						INT_4 maxLastNode, INT_4 minLastNode, INT_4 *aaPresentMass, INT_4 *seqNum, 
						INT_4 *subseqNum, SCHAR *sequenceNode);
int 			ExtensionsSortDescend(const void *n1, const void *n2);
struct extension	*SortExtensions(struct extension *inExtensionList);
struct extension 	Score2aaExtension(struct extension io2aaExtension, INT_4 startingNodeMass, 
						INT_4 *oneEdgeNodes, INT_4 oneEdgeNodesIndex,
						BOOLEAN oneEdgeNNode, char	sequenceNodeValue);
struct Sequence *AddExtensions(struct Sequence *subsequencePtr, SCHAR *sequenceNode, 
						INT_4 *oneEdgeNodes, INT_4 oneEdgeNodesIndex, 
						INT_4 *aaPresentMass, INT_4 topSeqNum, INT_4 *lastNode, INT_4 lastNodeNum, 
						INT_4 *seqNum, INT_4 maxLastNode, INT_4 minLastNode, INT_4 lowSuperNode, 
						INT_4 highSuperNode);
void 			FreeSequenceStructs(struct Sequence *s);
struct Sequence *AlterSubsequenceList(struct Sequence *firstPtr, struct Sequence *newPtr);
char 			CorrectMass(INT_4 *peptide, INT_4 peptideLength, INT_4 *aaPresentMass);
void 			amIHere(INT_4 correctPeptideLength, struct Sequence *subsequencePtr);
struct Sequence *LCQNterminalSubsequences(SCHAR *sequenceNode, INT_4 maxLastNode, INT_4 lowSuperNode,
						INT_4 highSuperNode);
void 			ClearLowestScore(struct Sequence *newSubsequencePtr, INT_4 *subseqNum, 
						INT_4 maxLastNode, INT_4 minLastNode);
				
/*Prototypes for LutefiskScore.*/
struct Sequence *ScoreSequences(struct Sequence *finalSequencePtr, 
						struct MSData *firstMassPtr);
void 			LoadTheIonArrays(struct MSData *firstMassPtr, INT_4 *fragNum, 
						INT_4 *fragMOverZ, INT_4 *fragIntensity);
INT_4 			TotalIntensity(INT_4 fragNum, INT_4 *fragMOverZ, 
						INT_4 *fragIntensity);
void 			WaterLoss(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, INT_4 *fragIntensity);
void 			LoadSequence(INT_4 *sequence, INT_4 *seqLength, 
						struct Sequence *currSeqPtr);
void 			PEFragments(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, 
						INT_4 *sequence, INT_4 seqLength);
char 			ArgIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, 
						INT_4 *sequence, INT_4 seqLength);
void 			InternalFrag(REAL_4 *ionFound, INT_4 *fragMOverZ, 
						INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *ionType);
INT_4 			FindABYIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				     	INT_4 *sequence, INT_4 seqLength, char argPresent,
				     	REAL_4 *yFound, REAL_4 *bFound, REAL_8 *byError, INT_4 *ionType);
void 			ScoreLowMassIons(REAL_4 *ionFound, INT_4 *fragMOverZ, 
						INT_4 *sequence, INT_4 seqLength, INT_4 lowMassIons[][3], INT_4 *ionType);
REAL_4 			IntensityScorer(INT_4 *fragIntensity, REAL_4 *ionFound, INT_4 cleavageSites, 
						INT_4 fragNum, INT_4 seqLength, INT_4 intensityTotal);
void 			FreeAllSequence(struct Sequence *currPtr);
void 			SeqIntensityRanker(struct SequenceScore *firstScorePtr);
void 			SeqComboScoreRanker(struct SequenceScore *firstScorePtr);
struct Sequence *AddTagBack(struct Sequence *firstSequencePtr);
void 			PrintToConsole(struct SequenceScore *firstScorePtr);
REAL_4			CalcIonFound(REAL_4 currentIonFound, INT_4 massDiff);
void			PrintToConsoleAndFile(struct SequenceScore *firstScorePtr, REAL_4 quality, INT_4 length, 
						REAL_4 perfectProbScore);
INT_4 			FindBYIons(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				     	INT_4 *sequence, INT_4 seqLength);
REAL_4 			BYIntensityScorer(INT_4 *fragIntensity, REAL_4 *ionFound, INT_4 cleavageSites, 
						INT_4 fragNum, INT_4 seqLength, INT_4 intensityTotal);
REAL_4 			IntensityOnlyScorer(INT_4 *fragIntensity, REAL_4 *ionFound, 
						INT_4 fragNum, INT_4 intensityTotal);
void 			CTerminalLysOrArg(struct MSData *firstMassPtr);
struct Sequence *InitLutefiskScore(INT_4 *sequence, INT_4 *fragMOverZ, INT_4 *fragIntensity,
						 REAL_4 *ionFound, REAL_4 *ionFoundTemplate, INT_4 *countTheSeqs,
						 struct Sequence *firstSequencePtr,
						 REAL_4 *yFound, REAL_4 *bFound, struct MSData *firstMassPtr,
						 INT_4 *charSequence, REAL_8 *byError, INT_4 *ionType);
void 			TossTheLosers(struct Sequence *firstSequencePtr, REAL_4 *ionFoundTemplate, INT_4 fragNum,
						INT_4 *fragMOverZ, INT_4 *fragIntensity, INT_4 intensityTotal,
						REAL_4 *ionFound, INT_4 *countTheSeqs, INT_4 *sequence);
INT_4 			FindNCharge(INT_4 *sequence, INT_4 seqLength);
char 			TwoAAExtFinder(INT_4 *sequence, INT_4 i);
INT_4 			BCalculator(INT_4 i, INT_4 *sequence, INT_4 bCalStart, INT_4 bCalCorrection);
INT_4 			YCalculator(INT_4 i, INT_4 *sequence, INT_4 seqLength, INT_4 yCalStart, 
						INT_4 yCalCorrection);
BOOLEAN 		CheckItOut(struct Sequence *firstSequencePtr);
INT_4 			AlterIonFound(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ,
				    	 INT_4 *sequence, INT_4 seqLength, REAL_4 *yFound, 
				    	 REAL_4 *bFound, INT_4 newCleavageSites);
void 			FreeAllSequenceScore(struct SequenceScore *currPtr);
void 			ScoreBYIsotopes(REAL_4 *ionFound, INT_4 *fragMOverZ, INT_4 fragNum, INT_4 *ionType);
void 			AdjustIonIntensity(INT_4 fragNum, INT_4 *fragIntensity);
void 			RevertBackToReals(struct MSData *firstMassPtr, struct SequenceScore *firstScorePtr);
BOOLEAN 		CheckItOutSequenceScore(struct SequenceScore *firstSequencePtr);
void 			HighMOverZFilter(struct Sequence *firstSequencePtr, INT_4 *fragMOverZ, 
						INT_4 *fragIntensity, INT_4 *countTheSeqs, 
						INT_4 *sequence, INT_4 fragNum);
REAL_4 			AssignHighMZScore(INT_4 highMZNum, INT_4 *highMZFrags, INT_4 *highMZInts, 
						REAL_4 totalIntensity, INT_4 *sequence, INT_4 seqLength);
struct Sequence *RemoveRedundantSequences(struct Sequence *firstSequencePtr);
char 			SingleAA(INT_4 aminoAcidMass);
char 			*PeptideString(INT_4 *peptide, INT_4 peptideLength);
REAL_4 			AssignError(REAL_4 currentError, INT_4 calculatedMass, INT_4 observedMass);
REAL_4 			StandardDeviationOfTheBYErrors(REAL_8 *byError, INT_4 fragNum);
void 			BoostTheCTerminals(struct MSData *firstMassPtr);
REAL_4 			ProInThirdPosition(REAL_4 oldScore, INT_4 *sequence, INT_4 seqLength);
INT_4			SequenceLengthCalc(INT_4 *sequence, INT_4 seqLength);
void 			MakeNewgGapList(void);
void			ExpandSequences(struct Sequence *firstSequencePtr);
char 			Ratchet(INT_4 *aaNum, INT_4 cycle, INT_4 *sequence, INT_4 seqLength);
char			GoodSequence(struct Sequence *firstSequencePtr, INT_4 *peptide, INT_4 *peptideCorrection, 
						INT_4 peptideLength);
void			AddToGapList(struct Sequence *firstSequencePtr);
char 			IsThisADuplicate(struct SequenceScore *firstScorePtr, INT_4 *sequence, 
						REAL_4 intOnlyScore, REAL_4 intScore, INT_4 seqLength);
REAL_4			Recalibrate(INT_4 fragNum, INT_4 *fragMOverZ, INT_4 *sequence, INT_4 seqLength,
						INT_4 *fragIntensity);
void 			ProlineInternalFrag(REAL_4 *ionFound, INT_4 *fragMOverZ, 
						INT_4 *sequence, INT_4 seqLength, INT_4 fragNum);
void			RescoreAndPrune(struct Sequence *firstSequencePtr, REAL_4 *ionFound, INT_4 fragNum, 
						INT_4 *fragMOverZ, INT_4 *sequence, INT_4 seqLength, 
						char argPresent, REAL_4 *yFound, REAL_4 *bFound, 
						REAL_8 *byError, INT_4 cleavageSites, 
						INT_4 lowMassIons[][3], REAL_4 *ionFoundTemplate,
						INT_4 *fragIntensity, INT_4 intensityTotal, INT_4 *ionType);
REAL_4			ScoreAttenuationFromCalfactor(REAL_4 calFactor, REAL_4 intScore);
void 			AddDatabaseSequences(struct Sequence *firstSequencePtr);
char 			*GetDatabaseSeq(INT_4 *peptide, INT_4 peptideLength);
void 			AddPyroGlu();
INT_4			FindNOxMet(INT_4 *sequence, INT_4 seqLength);
void 			RevertTheRevertBackToReals(struct MSData *firstMassPtr);
INT_4			SequenceLengthCalcNoFudge(INT_4 *sequence, INT_4 seqLength);
REAL_4			MassBasedQuality(INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *fragMOverZ, char argPresent);
BOOLEAN			KeepSequence(struct SequenceScore *currPtr, REAL_4 intscrKeep, REAL_4 xcorrKeep, 
						REAL_4 qualityKeep, REAL_4 probscrKeep);
REAL_4			ComboScore(struct SequenceScore *currPtr);

struct SequenceScore *LoadSeqScoreStruct(REAL_4 intScore, REAL_4 intOnlyScore,
						INT_4 *sequence, INT_4 *charSequence, INT_4 seqLength,
						REAL_4 stDevErr, INT_4 cleavageSites, REAL_4 calFactor,
						char databaseSeq, REAL_4 normalXCorScore, REAL_4 quality,
						REAL_4 length, REAL_4 probScore, REAL_4 comboScore);
struct SequenceScore *AddToSeqScoreList(struct SequenceScore *firstPtr, 
											struct SequenceScore *currPtr);
struct SequenceScore *FindLowestScore(struct SequenceScore *currPtr);
struct SequenceScore *MassageScores(struct SequenceScore *firstScorePtr);
INT_4 ScoreC1(REAL_4 *ionFound, INT_4 fragNum, INT_4 *fragMOverZ, 
				 INT_4 *sequence, INT_4 seqLength);
						
/*Prototypes for LutefiskXCorr.*/
extern void 	DoCrossCorrelationScoring(struct SequenceScore *firstScorePtr, struct MSData *firstMassPtr) ;
void 			CrossCorrelate(REAL_4 *array1, REAL_4 *array2, UINT_4 n, REAL_4 *result);
REAL_4 			YXCorrCalc(INT_4 i, struct SequenceScore *currScorePtr, REAL_4 YionStart, 
						INT_4 seqLength);
REAL_4 			BXCorrCalc(INT_4 i, struct SequenceScore *currScorePtr, REAL_4 BionStart);
INT_4 			FindNChargeXCorr(struct SequenceScore *currScorePtr);
void 			AddPeakToSpectrum( REAL_4 *spectrum, REAL_4 mass, REAL_4 intensity);

/*Prototypes for LutefiskFourier.*/                                      
void 			CalcNormalizedExptPeaks(struct MSData *firstMassPtr);
void 			FillInSpectrum1(struct MSData *firstMassPtr);
void 			SetupCrossCorrelation(void);

/*Prototypes for LutefiskGetAutoTag.*/
void 			MaskSequenceNodeWithTags(SCHAR *sequenceNode, char *tagNode);
INT_4 			TagExtensions(char *tagNode, INT_4 topSeqNum, INT_4 *tagNodeIntensity,
						INT_4 totalIonIntensity);
INT_4 			NterminalTags(char *tagNode, INT_4 *tagNodeIntensity);
void 			TagMaker(char *tagNode, INT_4 *tagNodeIntensity,
						INT_4 totalIonIntensity);
void 			FindTagYIons(char *tagNode, INT_4 charge, INT_4 *tagNodeIntensity);
void 			TagNodeInit(char *tagNode, INT_4 *tagNodeIntensity);
void 			GetAutoTag(struct MSData *firstMassPtr, SCHAR *sequenceNode);
struct Sequence *AlterTagList(struct Sequence *firstPtr, struct Sequence *newPtr);
void 			FreeTagStructs(struct Sequence *currPtr);
struct Sequence *LoadFinalTagStruct(INT_4 *peptide, INT_4 peptideLength, 
						INT_4 score, INT_4 nodeValue, INT_4 gapNum);
struct Sequence *LinkTagList(struct Sequence *firstPtr, struct Sequence *newPtr);
struct Sequence *LoadTagStruct(INT_4 *peptide, INT_4 peptideLength, 
						INT_4 score, INT_4 nodeValue, INT_4 gapNum);
void 			FreeMSData(struct MSData *currPtr);
struct MSData 	*LoadMSData(REAL_4 massValue, INT_4 ionIntensity);
struct MSData 	*AddCIDList(struct MSData *firstPtr, struct MSData *currPtr);
void 			FindTagB2Ions(struct MSData *firstMassPtr, INT_4 *totalIntensity, 
						char *tagNode, INT_4 *tagNodeIntensity);
void 			FindMoreB2Ions(struct MSData *firstMassPtr, INT_4 *totalIntensity, 
						char *tagNode, INT_4 *tagNodeIntensity);
INT_4 			AlternateNterminalTags(char *tagNode, INT_4 *tagNodeIntensity);
void 			AlternateTagMaker(char *tagNode, INT_4 *tagNodeIntensity,
						INT_4 totalIonIntensity);
void 			FilterTagMasses(INT_4 charge);
int 			SequenceScoreDescendSortFunc(const void *n1, const void *n2);
char 			*ComposePeptideString(INT_4 *peptide, INT_4 peptideLength);
void 			RemoveNeutralLosses(INT_4 charge);
void 			PrintScoreDetailsToXLFile(struct SequenceScore *firstScorePtr, REAL_4 perfectProbScore);
REAL_4			CalcPerfectProbScore(INT_4 fragNum, INT_4 *fragMOverZ);
struct SequenceScore 	*DetermineBestCandidates(struct SequenceScore *firstScorePtr);
struct extension 		*SortExtension(struct extension *inExtensionList);

/*Prototypes for LutefiskProbScorer*/
REAL_4			LutefiskProbScorer(INT_4 *sequence, INT_4 seqLength, INT_4 fragNum, INT_4 *fragMOverZ, char argPresent);
REAL_8			FindImmoniumIons(INT_4 *mass, INT_4 ionCount, 
						REAL_8 probScore, REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength);
REAL_8			FindInternalIons(INT_4 *mass, INT_4 ionCount, 
						REAL_8 probScore, REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength);
REAL_4			InitProbScore(INT_4 *sequence, INT_4 seqLength);
REAL_8 			FindBIons(INT_4 *mass, INT_4 ionCount, REAL_8 probScore,
						REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength, char argPresent);
REAL_8 			FindYIons(INT_4 *mass, INT_4 ionCount, REAL_8 probScore,
						REAL_4 *randomProb, INT_4 *sequence, INT_4 seqLength, char argPresent);
void	 		CalcRandomProb(REAL_4 *randomProb, INT_4 *mass, INT_4 ionCount);


#endif /*	__LUTEFISK__*/


