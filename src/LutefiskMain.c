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
*********************************************************************************************
        
   Lutefisk is a program designed to aid in the interpretation of CID data of peptides.  
   The main assumptions are that the data is of reasonable quality, the N- and C-terminal
   modifications (if any) are known, and the precursor ion charge (and therefore the 
   peptide molecular weight) are known.  The ultimate goal here is to develop code that
   can utilize msms data in conjunction with ambiguous and incomplete Edman sequencing data,
   sequence tags, peptide derivatization, and protein or est database searches.  An older 
   version of Lutefisk has been written in FORTRAN and runs on 68K Macs that have an fpu
   (1991, 39th ASMS Conference on Mass Spectrometry and Allied Topics, Nashville, TN, pp 1233-
   1234).  This is a different and improved algorithm partly inspired by Fernandez-de-Cossjo, 
   et al. (1995) CABIOS Vol. 11 No. 4 pp 427-434.  Combining this msms interpretation algorithm
   with Edman sequencing, database searches, and derivatization is entirely of my own design;
   J. Alex Taylor implemented the changes in the FASTA code (Bill Pearson, U. of VA) so that
   the Lutefisk output can be read directly by the modified FASTA program.  In addition, 
   there were a number of additional critical changes made to FASTA to make it more compatible 
   with msms sequencing data.
        
   The trademark Lutefisk was chosen at random, and is not meant to imply any similarity
   between this computer program and the partially base-hydrolzyed cod fish of the same name. 
********************************************************************************************/

/* ANSI headers */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#if(defined(__MWERKS__) && __dest_os == __mac_os)
/* Some Macintosh specific things */
    #include "getopt.h"

    #include <StandardFile.h>		
    #include <sioux.h>
    #include <console.h>			
StandardFileReply freply;   
Point wpos;         
INT_4 tval;         
char prompt[256];       
#endif

#if(defined(__MWERKS__) && __dest_os == __win32_os)

/* Some Windoze specific things */
    #include "getopt.h"
#endif


/* Lutefisk headers */
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"


char versionString[256] = "LutefiskXP v1.0.7\nCopyright 1996-1906 Richard S. Johnson\n\n";       

/*
//--------------------------------------------------------------------------------
//  main()
//--------------------------------------------------------------------------------
  There are six main steps to be performed by Lutefisk:

  1. Import the information from the editable text file Lutefisk.params and Lutefisk.details.
     Lutefisk.params contains information like the CID data file name, the peptide molecular
     weight, and so on.  Lutefisk.details contains more esoteric information that will
     not typically be modified, but that would be nice to be able to alter.  For example,
     the .details file contains info relating to the types of ions to be considered, and the 
     scoring values for those ions.  In addition, there is the file Lutefisk.edmans, 
     which contains the Edman data, if its available.
  2. Import the CID data.  This is currently done by making an ASCI file of m/z versus
     intensity using Finnigan's "List" program.  Then this ASCI file is imported by
     Lutefisk.  Its also possible to import a tab-delimited ASCI file not produced by Finnigan.
  3. Find the nodes on the sequence graph.  The CID data is converted to a sequence graph, 
     which is where each ion is assumed to be one of the ion types under consideration
     (b, or y, or whatever) and then mathematically converted to the corresponding b ion
     value.  For example, an ion at m/z 800 for a peptide of molecular weight 999 when
     assumed to be a y ion would mathematically be equivalent to a b ion of m/z 200.
     The nodes on this sequence graph are integer values, so that differences between
     nodes are directly comparable to nominal mass values of amino acid residues.  Each
     node will have two scores - Nterm and Cterm - derived from the number of ions that
     provide evidence for each node (C-terminal ions and N-terminal ions).  If a sequence 
     tag is available, then it can be used to reduce the size of the sequence graph in that 
     those nodes corresponding to the tag are removed.  There are three ways to make the sequence 
     graph: a general peptide graph, a tryptic triple quad graph, and a tryptic ion trap graph.
  4. Calculate the Summed Node Scores starting from the C-terminus, and identify all of the
     one-edged nodes.  Beginning with the C-terminus, start connecting the nodes, and assign to
     each connected node an additional bonus score.  One edge nodes are those that are connected
     to the C-terminus but cannot be extended toward the N-terminus.  These one-edge nodes 
     are used in step #5 for jumping two amino acid gaps in the fragmentation pattern.
  5. Build the subsequences beginning at the N-terminus, and keep those with the highest
     subsequence scores that are derived from the node values established in step #4.  
     If a subsequence cannot be extended, then the program looks for any two amino acid jumps 
     to one of the one-edge nodes from step #4.  Unmodified peptides begin at mass 1 (hydrogen), 
     otherwise they start at 43 (acetylation), 44 (carbamylation), or 112 (pyroglutamic acid).
     In any case, the first jump from the N-terminus can involve either one or two amino acids. 
     Thereafter, only one amino acid is added at a time.  For ion trap data, the initial jump can 
     be more than two amino acids.
  6. Score the completed sequences obtained from step #5.  The score will be the usual
     percentage of ion current accounted for, plus it might be good to add the option
     of including the cross-correlation score.

*************************************************************************************************/

int main(int argc, char **argv)
{
    INT_4 i;
    const   time_t          theTime = (const time_t)time(NULL);
    extern INT_4 optind;

#if(defined(__MWERKS__) && __dest_os == __mac_os)
    argc = ccommand(&argv);  
#endif

/*    gParam.startTicks = clock();*/

    if (!SystemCheck()) exit(1);

    BuildPgmState(argc, argv);      /*Read in any line commands*/

    if (gParam.fMonitor)
    {
        printf(versionString);
        printf("Run Date: %20s", ctime(&theTime));
    }

    ReadParamsFile();


    ReadResidueFile();

/*  
*       Import the information from Lutefisk.details.
*
*       The various ion values are read from Lutefisk.details and are used to assign values to the
*       nodes in the sequence graph.  The value "fragmentPattern" is needed in order to figger out
*       which of the three columns should be read in as the set of ion values.  The function returns
*       a value that corresponds to the sum of all of the ion values for b, y, etc ions.
*/

    gIonTypeWeightingTotal = ReadDetailsFile();

/*      Reassign values for Cys to account for changes in mass due to alkylation.*/

    gMonoMass[C] = gParam.cysMW;
    gAvMass[C]   = gParam.cysMW;
    gNomMass[C]  = gParam.cysMW;



/*      Assign values to the globals H2O and NH3.*/
    H2O = 2 * gElementMass[HYDROGEN] + gElementMass[OXYGEN];
    NH3 = gElementMass[NITROGEN] + 3 * gElementMass[HYDROGEN];    

/*
*	The Edman data is read from the file Lutefisk.edman if the parameter "edmanPresent" equals
*	'Y'.  The data is read into an array gEdmanData[cycle number][amino acids in the cycle], and
*	it contains the nominal mass values of the amino acids found in the file Lutefisk.edman
*	(rather than a character listing of the amino acid single letter code).
*/

    if (gParam.edmanPresent)
    {
        ReadEdmanFile();
    }


    /* Total hack because we mess with these values later on. */
    gParam.peptideMW_orig   = gParam.peptideMW;
    gParam.peptideErr_orig  = gParam.peptideErr;
    gParam.fragmentErr_orig = gParam.fragmentErr;
    gParam.peakWidth_orig   = gParam.peakWidth;
    gParam.monoToAv_orig    = gParam.monoToAv;
    gParam.qtofErr_orig     = gParam.qtofErr;
    gParam.ionOffset_orig   = gParam.ionOffset;
    gParam.cysMW_orig       = gParam.cysMW;
    gParam.tagNMass_orig    = gParam.tagNMass;
    gParam.tagCMass_orig    = gParam.tagCMass; 
    gParam.maxGapNum_orig   = gParam.maxGapNum; 
    gParam.modifiedNTerm_orig = gParam.modifiedNTerm;
    gParam.modifiedCTerm_orig = gParam.modifiedCTerm;
    gParam.topSeqNum_orig	= gParam.topSeqNum;
    strcpy(gParam.outputFile_orig, gParam.outputFile);

//	optind = 0;	/*debug*/
//	argc = 3;	/*debug*/
    if (optind < argc)
    {
        for (i = optind; i < argc; i++)
        {
            strcpy(gParam.cidFilename, argv[i]);
            Run();
        }
    }
    else if (strlen(gParam.cidFilename) > 0)
    {
        /* A filename was specified in the Lutefisk.params file */
       
        Run();
        
    }

    sleep(1);   /* Why sleep before quitting? Well, on our 500 MHz alpha some results 
                   were being lost when calling the program from a child process via a pipe; 
                   seemingly because the pipe was being terminated before all the data had
                   gotten through. */
    return(0);        /* All done */

}
/*************************************************************************************************/
void Run()
{
    REAL_4 actualPeptideMW, actualTopSeqNum, actualFinalSeqNum; 
    INT_4 i, massChange, posNeg;
    SCHAR *sequenceNode = NULL;
    SCHAR *sequenceNodeC = NULL;
    SCHAR *sequenceNodeN = NULL;
    INT_4 oneEdgeNodesIndex, *oneEdgeNodes = NULL;

    struct MSData *firstMassPtr = NULL, *firstRawDataPtr = NULL;
    struct Sequence *firstSequencePtr = NULL;
    const   time_t          theTime = (const time_t)time(NULL);

	gParam.startTicks = clock();
    /* Total hack because we mess with these values later on. */
    gParam.peptideMW   = gParam.peptideMW_orig;
    gParam.peptideErr  = gParam.peptideErr_orig;
    gParam.fragmentErr = gParam.fragmentErr_orig;
    gParam.peakWidth   = gParam.peakWidth_orig;
    gParam.monoToAv    = gParam.monoToAv_orig;
    gParam.qtofErr     = gParam.qtofErr_orig;
    gParam.ionOffset   = gParam.ionOffset_orig;
    gParam.cysMW       = gParam.cysMW_orig;
    gParam.tagNMass    = gParam.tagNMass_orig;
    gParam.tagCMass    = gParam.tagCMass_orig; 
    gParam.maxGapNum   = gParam.maxGapNum_orig; 
    gParam.modifiedNTerm = gParam.modifiedNTerm_orig;
    gParam.modifiedCTerm = gParam.modifiedCTerm_orig;
    gParam.topSeqNum	= gParam.topSeqNum_orig;
    strcpy(gParam.outputFile, gParam.outputFile_orig);

    gFirstTimeThru = TRUE;
/*
*       GetCidData opens an ASCII file containing lists of m/z values and intensities for the 
*       CID data.  This file is produced using the Finnigan program called "LIST", by using the 
*       "Print..." command found in the "File" menu.  A dialog box appears where you tell it to
*       the saved format is "ASCII" for both the text and graph displays, and then you provide
*       a file name and click on the "save to file" button.
*       The pointer firstMassPtr points to the first
*       element in the linked list of ion values (m/z and intensity).
*
*       GetCidData uses gElementMass and gMonoMass instead of gElementMass_x100 and gMonoMass_x100,
*       which haven't even been assigned values yet.
*/      

    firstMassPtr = GetCidData();
    if (NULL == firstMassPtr)
    {
        PrintPartingGiftToFile();
        exit(0);
    }    

/*
*	Adjust peptideMW for LCQ data using fragment ion pairs.  
*/   

	if(gParam.fragmentPattern == 'L')
	{
		AdjustPeptideMW(firstMassPtr);
	}

/*
*       Change output filename if no output filename specified - start name plus ".lut".
*       martin 98/8/27
*/

    ChangeOutputName();

/*
*       If the peptideMW is obtained from the data file header, then the sequence tag cannot be
*       set up properly until peptideMW is known.  Now that the data file has been read, I can
*       set this up correctly.
*/

    SetupSequenceTag();


/*      Assign space to the various arrays.*/

    sequenceNode = (SCHAR *) malloc(gGraphLength * sizeof(SCHAR ));    /*Will contain summary of evidence.*/
    if (sequenceNode == NULL)
    {
        printf("main:  Out of memory");
        exit(1);
    }
    sequenceNodeC = (SCHAR *) malloc(gGraphLength * sizeof(char ));   /*Will contain C-terminal evidence.*/
    if (sequenceNodeC == NULL)
    {
        printf("main:  Out of memory");
        exit(1);
    }
    sequenceNodeN = (SCHAR *) malloc(gGraphLength * sizeof(char ));   /*Will contain N-terminal evidence.*/
    if (sequenceNodeN == NULL)
    {
        printf("main:  Out of memory");
        exit(1);
    }
    oneEdgeNodes = (int *) malloc(gGraphLength * sizeof(INT_4 ));   /*Will contain evidence that only
                                                                      connects w/ the C-terminus.*/
    if (oneEdgeNodes == NULL)
    {
        printf("main:  Out of memory");
        exit(1);
    }


/*      If the gParam.maxGapNum is equal to -1, then assign a value based on gParam.peptideMW.*/

    if (gParam.maxGapNum == -1)
    {
        if (gParam.peptideMW < 1400)
        {
            gParam.maxGapNum = 1;
        }
        else if (gParam.peptideMW >= 1400 && gParam.peptideMW < 2000)
        {
            gParam.maxGapNum = 2;
        }
        else if (gParam.peptideMW >= 2000)
        {
            gParam.maxGapNum = 3;
        }
    }

/*
*       Multiply the gElementMass and gMonoMass values to give integer numbers for the corresponding
*       arrays of gElementMass_x100 and gMonoMass_x100.  These latter arrays are used to represent the 
*       fractional mass values of the elemental and amino acid masses.  The defined value of GRAPH_LENGTH
*       is divided by 10, 100, 1000, 10,000 and 100,000 until a value less than 10,000 is obtained.  Also,
*       this is where gParam fields that are mass-related are multiplied by the gMultiplier value.
*/

    CreateGlobalIntegerMassArrays(firstMassPtr);

/*
*       Make and modify the array gGapList so that it incorporates information about amino acids
*       that are missing, plus it alters the mass of cysteine based on the value of cysMW.  The first
*       gAminoAcidNumber positions in gGapList contain the nominal residue mass values for the amino
*       acids, except that the masses of Gln and Ile are assigned zero.  Hence, I need to be careful
*       about this throughout the program, and I often have an if statement that prevents using
*       gGapList values of zero.  Positions after gAminoAcidNumber contain residue masses for 
*       two amino acids.
*
*       The other global arrays (gNomMass, gSingAA, etc.) remain intact.  Later, in
*       ScoreSequences I modify some globals that are global only to the functions within that file.
*
*       SetupGapList uses gElementMass_x100 and gMonoMass_x100.
*/

    SetupGapList();


/*
*       MakeSequenceGraph assigns values to the array of INT_4 's sequenceNodeC[GRAPH_LENGTH]
*       and sequenceNodeN[GRAPH_LENGTH].
*       The indexing of these arrays corresponds to nominal mass values of hypothetical b-type
*       ions, and the values assigned to each node is an estimation of the likelihood that 
*       there is a real cleavage at that mass.
*
*       MakeSequenceGraph uses gElementMass_x100 and gMonoMass_x100.
*
*       Here's where I would start looping at different pParam.peptideMW values in order to 
*       determine average scores for incorrect sequences.  
*/
    for (i = 0; i <= gParam.wrongSeqNum; i++)    /*initialize*/
    {
        gWrongXCorrScore[i]		= 0;	/*holds the best wrong cross-correlation scores*/
        gWrongIntScore[i]		= 0;	/*holds the best wrong intensity scores*/
        gWrongProbScore[i]		= 0;	/*holds the best wrong Pevzner score*/
        gWrongQualityScore[i]	= 0;	/*holds the best wrong quality*/
        gWrongComboScore[i]		= 0;	/*holds the best combined score*/
    }

    actualPeptideMW = gParam.peptideMW;     /*save the real peptide mass*/
    actualTopSeqNum = gParam.topSeqNum;     /*save the real max subsequence number*/
    actualFinalSeqNum = gParam.finalSeqNum; /*save the real max final candidate sequence number*/
    posNeg = 1; /*goes back and forth between +1 and -1, see below*/

    /*Here's the loop where i is negative and works towards zero, which represents the correct
      mass.  If i is -10 then -9, massChange becomes -5 and -5.  However, each time thru the 
      loop posNeg goes back and forth from +1 to -1, so in the end massChange is -5 then +5.  
      These are the numbers that get multiplied by a methylene mass and then added to the correct
      peptide mass.*/
    for (i = -1 * gParam.wrongSeqNum; i <= 0; i++)
    {
        if (i != 0)
        {
            massChange = (REAL_4)i / 2 - 0.5;   /*since i is neg I need to subtract 0.5 to round
                                                    down*/
            massChange = massChange * posNeg;   /*go pos and neg*/
            posNeg = posNeg * -1;   /*go back and forth -1 +1 -1 +1 on and on*/
            gCorrectMass = FALSE;
            /*for wrong masses use a smaller seqNum to speed the processing*/
            gParam.topSeqNum = 1000;
            gParam.finalSeqNum = 5000;
        }
        else
        {
            massChange = 0; /*no mass change when i = 0, cuz thats the loop for the correct MW*/
            gCorrectMass = TRUE;    /*loop is for correct mass*/
            gParam.topSeqNum = actualTopSeqNum; /*use correct subseq num for correct mass*/
            gParam.finalSeqNum = actualFinalSeqNum;/*use correct final sequence num */
        }

        /*gParam.peptideMW gets changed for the remainder of the loop*/
        gParam.peptideMW = actualPeptideMW + massChange * 
                           (/*2 * gElementMass_x100[HYDROGEN]*/ + gElementMass_x100[CARBON]);
							/*differences of a methylene is debatable; I think it might be bad idea now*/

        MakeSequenceGraph(firstMassPtr, sequenceNode, sequenceNodeC, sequenceNodeN, 
                          gIonTypeWeightingTotal);

/*
*       SummedNodeScore connects the nodes starting from the C-terminal node(s) that differ by the
*       nominal mass of an amino acid residue.  There may be several C-terminal nodes if the peptide
*       mass error is sufficiently large, all of which are used independently of each other.  Those
*       nodes that can be connected to the C-terminus are given a bonus score to differentiate them
*       from those nodes that do not.   This is also
*       where the C-terminal one-edge nodes are found and stored in the array oneEdgeNodes, and
*       has a maximum of oneEdgeNodesIndex (ie, if oneEdgeNodesIndex is 26, then oneEdgeNodes has values
*       for [0 to 25]).  The altered node values are held in the array sequenceNode.  The arrays
*       sequenceNodeN and sequenceNodeC were obtained from the function MakeSequenceGraph and are used
*       as the input information for SummedNodeScore.  At one time, this function summed the node
*       scores that lead up to it, but I found that this tended to overly dominate the subsequencing
*       scores in a bad way.  What worked best is to add a bonus score to any node that connects to the
*       C-terminus.  However, the name of the function remains - SummedNodeScore - even though it
*       doesn't sum the node scores.  Those bits of evidence in either SequenceNodeC or SequenceNodeN
*       that cannot connect to the C-terminus are added together to give a relatively low node
*       value.  Also, this is only done for stretches of consecutive nodes where none of them
*       are able to connect to the C-terminus (ie, nodes adjacent to one that does connect are not
*       included in the final array sequenceNode).
*
*       SummedNodeScore uses gElementMass_x100 and gMonoMass_x100.
*/

        SummedNodeScore(sequenceNode, sequenceNodeC, sequenceNodeN, oneEdgeNodes,
                        &oneEdgeNodesIndex, gIonTypeWeightingTotal);


/*
*       GetAutoTag finds bits of sequences using only the m/z region between the precursor ion and
*       above.  Ions of type y are only considered unless there are pairs of ions that differ by 28
*       where the higher m/z pair is of greater intensity.  The only charge state considered is
*       one less than the precursor charge.  For LCQ data, where both b and y ions are usually seen
*       above the precursor, this procedure can be a bit dangerous.  I've found situations where
*       the correct sequence is eliminated because the correct sequence tag is not found.  Due
*       to the ambiguity of LCQ data w/ respect to it being a b or a y, I don't use the auto-tag
*       feature for trap data.  It works great for TSQ data, though.
*
*       GetAutoTag uses gElementMass_x100 and gMonoMass_x100.
*/

        if ((gParam.fragmentPattern == 'L' || gParam.fragmentPattern == 'T' || gParam.fragmentPattern == 'Q')
            && gParam.chargeState > 1 && gParam.autoTag)
        {
            GetAutoTag(firstMassPtr, sequenceNode);
        }


/*      
*       Now that the node scores have been finalized (in the array sequenceNode) and the C-terminal
*       one-edged nodes have been identified, it is time to start building up subsequences from the
*       N-terminus.  Again, I connect the nodes that are spaced one or two amino acid residues apart,
*       but now I need to remember how the nodes were connected.  In SummedNodeScore, all that was 
*       important was knowing that there was some way to connect the nodes from the C-terminus, whereas
*       here I need to keep track of the pathway from the N-terminus.  The output is the final list of
*       completed sequences, which is contained in a struct of type Sequence.  The function 
*       SubsequenceMaker returns a pointer to a struct of type Sequence, which is the first element in
*       a linked list of completed sequences plus the associated subsequence score.  If there were
*       no completed sequences, then the function returns a NULL value.
*
*       SubsequenceMaker uses gElementMass_x100 and gMonoMass_x100.
*/

        firstSequencePtr = SubsequenceMaker(oneEdgeNodes, oneEdgeNodesIndex, sequenceNode);
        
/*
*		Next add subsequences that do not necessarily connect to either termini (process called Haggis).  
*		Mass scrambles for statistics has to be turned off, since changes in peptide MW will not alter the 
*		results.  
*/
		
		if(gParam.wrongSeqNum == 0)
		{
			firstSequencePtr = Haggis(firstSequencePtr, firstMassPtr);
		}

/*
*       Next, the list of sequences in the linked list starting w/ firstSequencePtr are scored. 
*       To do this, I need the CID data (firstMassPtr), the peptide molecular weight, the fragment
*       ion error, the charge state of the precursor ion, and the mass of cysteine.  The latter is
*       used to determine if certain alkylating groups are present that give rise to specific types
*       of fragment ions.  The return is a pointer to the ranked and scored list of sequences.
*       The linked list starting with firstMassPtr is intact, but all lists of sequences, except
*       for the returned linked list, are free'ed.  In addition, I need to know the monoToAv mass
*       switch, the sequence tag information (to add the sequence tag back in), and the N-terminal
*       modification.  Intensity-based scores are determined for each sequence, and the the top
*       MAX_X_CORR_NUM sequences are assigned cross-correlation scores.  In the end, some combination
*       of these two scores will provide a INT_4 list of sequences to be submitted for FASTA database
*       analysis.
*
*       ScoreSequences uses gElementMass and gMonoMass instead of the _x100 arrays.
*/


        if (firstSequencePtr != NULL)
        {
            firstSequencePtr = ScoreSequences(firstSequencePtr, firstMassPtr);
            gFirstTimeThru = FALSE; /*forever false after first time thru loop*/
            SetupGapList(); /*The gGapList can get changed in the scoring, so its returned to the 
                            original values*/
            FreeSequence(firstSequencePtr); /*Get rid of the sequences, cuz you'll be getting
                                            a new set soon*/
        }
        else if (gCorrectMass)
        {
            PrintPartingGiftToFile();
        }
        else
        {
            gWrongIndex++; /*No sequences to score, so the best score is zero, which is what
                     the gScore arrays were normalized to*/
        }


    }   /*end of gParam.peptideMW looping*/

    /*trash these things*/
    free(sequenceNodeC);    
    free(sequenceNodeN);
    free(oneEdgeNodes);
    free(sequenceNode);

/*      Free up the linked lists.*/
/* JAT - Why not free if win32? */
#if (__dest_os != __win32_os)

    /*List of completed sequences from subsequencing routine*/
    FreeMassList(firstMassPtr);             /*List of ions and intensities*/
#endif

    fflush(stdout);
}


/*
//--------------------------------------------------------------------------------
//  BuildPgmState()
//--------------------------------------------------------------------------------
*/

static void BuildPgmState(INT_4 argc, CHAR **argv)
{

    INT_4 c;

    extern CHAR *optarg;
    extern INT_4 optind;


    /* initialize parameters */

    gParam.fMonitor = TRUE;

    gParam.fVerbose = TRUE;

    strcpy(gParam.paramFile,"Lutefisk.params");

    strcpy(gParam.outputFile,"");

    strcpy(gParam.detailsFilename,"Lutefisk.details");

    strcpy(gParam.residuesFilename,"Lutefisk.residues");



    /* get command-line parameters */

    while ((c = getopt(argc, argv, "?hqvd:o:m:p:r:s:")) != -1)
    {

        switch (c)
        {
        
        case 'o':
            /* output file name */
            strncpy(gParam.outputFile, optarg, sizeof(gParam.outputFile));
            break;

        case 'd':
            /* details file name */
            strncpy(gParam.detailsFilename, optarg, sizeof(gParam.detailsFilename));
            break;  

        case 'm':
            /* peptide MW */
            gParam.peptideMW = atof(optarg);
            break;

        case 'p':
            /* param file name */
            strncpy(gParam.paramFile, optarg, sizeof(gParam.paramFile));
            break;

        case 'q':
            /* QUIET! */
            gParam.fMonitor = FALSE;
            gParam.fVerbose = FALSE;
            break;

        case 'r':
            /* residues file name */
            strncpy(gParam.residuesFilename, optarg, sizeof(gParam.residuesFilename));
            break;  

        case 's':
            /* database sequences file */
            strncpy(gParam.databaseSequences, optarg, sizeof(gParam.databaseSequences));
            break;

        case 'v':
            /* verbose */
            gParam.fVerbose = TRUE;
            break;

        case '?':
        case 'h':
            /* print usage */
            puts("\nUSAGE:  lutefisk [options] [CID file pathname]\n");
            puts(  "                -o = output file pathname");
            puts(  "                -q = quiet mode ON (default OFF)");
            puts(  "                -m = precursor ion mass");
            puts(  "                -d = details file pathname");
            puts(  "                -p = params file pathname");
            puts(  "                -r = residues file pathname");
            puts(  "                -s = pathnane of file with database sequences to score");
            puts(  "                -v = verbose mode ON (default OFF)");
            puts(  "                -h = print this help text");
            puts(  "" );
            puts("\n");
            exit(1);

            break;
        }
    }


    /* report flag state */

    if (gParam.fMonitor)
    {
        printf("Verbose mode %s\n",
               gParam.fVerbose ? "ON" : "OFF");
    }


    /* get the parameter file name if one is specified */

    /*XXXX       if (optind < argc) 
           {
                   strcpy(gParam.paramFile, argv[argc-1]);
           }
   */


    return;

}

/*
//--------------------------------------------------------------------------------
//  FindTheMultiplier()
//-------------------------------------------------------------------------------
    FindTheMultiplier uses the value of gParam.fragmentErr to determine the value of 
    gMultiplier. This is also where a value is determined for GRAPH_LENGTH, which will 
    replace the #defined value of GRAPH_LENGTH.
*/

void FindTheMultiplier(void)
{
    INT_4 multiplier;
    REAL_4 testMass;

    multiplier = 1;
    testMass = multiplier * gParam.fragmentErr;
    if (testMass > MULTIPLIER_SWITCH)        /* A value of two implies that there are a total of 5 nodes (2 on each side)*/
    {
        gMultiplier = 1;
    }
    else
    {
        multiplier = 10;
        testMass = multiplier * gParam.fragmentErr;
        if (testMass > MULTIPLIER_SWITCH)
        {
            gMultiplier = 10;
        }
        else
        {
            multiplier = 100;
            testMass = multiplier * gParam.fragmentErr;
            if (testMass > MULTIPLIER_SWITCH)
            {
                gMultiplier = 100;
            }
            else
            {
                multiplier = 100;
                testMass = multiplier * gParam.fragmentErr;
                if (testMass > MULTIPLIER_SWITCH)
                {
                    gMultiplier = 100;
                }
                else
                {
                    multiplier = 1000;
                    testMass = multiplier * gParam.fragmentErr;
                    if (testMass > MULTIPLIER_SWITCH)
                    {
                        gMultiplier = 1000;
                    }
                    else
                    {
                        gMultiplier = 1000;
                    }
                }
            }
        }
    }

/*      Now calculate GRAPH_LENGTH, which will a bit larger than required for the peptide mass.*/



    gGraphLength = gMultiplier * (gParam.peptideMW + gParam.wrongSeqNum 
                                  * (2 * gElementMass[HYDROGEN] + gElementMass[CARBON])) * 1.1;

    return;
}

/*
//--------------------------------------------------------------------------------
//  ChangeOutputName()
//--------------------------------------------------------------------------------
    If no output filename specified, set the output filename to the CID filename + ".lut".
    martin 98/8/27
    
    modified 000310 JAT
*/

void ChangeOutputName(void)
{
    if (strlen(gParam.outputFile) == 0)
    {
        /* There wasn't an output filename specified. */
        char  outputFile[256];
        INT_4 length;
        INT_4 fileCount;

        /* Start from the CID filename */
        strcpy (outputFile, gParam.cidFilename);

        length = strlen(outputFile);

        /* Add ".lut" to the end of the name (replacing .dta, etc.) */
        if ((length > 4) 
            && (!strncmp(outputFile + length - 4, ".dta", 4)
                || !strncmp(outputFile + length - 4, ".DTA", 4)
                || !strncmp(outputFile + length - 4, ".dat", 4)
                || !strncmp(outputFile + length - 4, ".DAT", 4)
                || !strncmp(outputFile + length - 4, ".txt", 4)
                || !strncmp(outputFile + length - 4, ".TXT", 4))
           )
        {
            strcpy(outputFile + length - 4, ".lut");
        }
        else
        {
            strcat(outputFile, ".lut");
        }


        /* Make sure that the file doesn't already exist. If it does, append a number. */
        strcpy(gParam.outputFile, outputFile);
        fileCount = 1;

        while (1)
        {
            FILE *fp = fopen(gParam.outputFile, "r");

            if (NULL == fp) break;

            fclose(fp);

            strcpy(gParam.outputFile, outputFile);
            sprintf(gParam.outputFile + strlen(gParam.outputFile), "%d\0", fileCount++);

            if (fileCount > 20)
            {
                printf("Too many old output files! Please clean up a bit first! Quitting.");
                exit(1);
            }
        }     


    }

    return;
}

/*
//--------------------------------------------------------------------------------
//  CreateGlobalIntegerMassArrays()
//--------------------------------------------------------------------------------
    The value of GRAPH_LENGTH is divided by 10, 100, 1000, 10000, and 100000 until 
    a value less than 10,000 is obtained.  The divisor obtained is then used to multiply
    the float values in gElementMass and gMonoMass (correctly rounded) to be placed in 
    the arrays gElementMass_x100 and gMonoMass_x100 (so-named because initially I will 
    use monoisotopic masses down to 2 decimal points).
*/
void CreateGlobalIntegerMassArrays(struct MSData *firstMassPtr)
{
    INT_4 i;
    REAL_4 correction;
    struct MSData *currPtr;

/*      Create integer values for the monoisotopic masses of amino acids.*/

    for (i = 0; i < gAminoAcidNumber; i++)
    {
        gMonoMass_x100[i] = (gMonoMass[i] * gMultiplier) + 0.5;
    }

/*      Create integer values for gNodeCorrection.*/

    for (i = 0; i < gAminoAcidNumber; i++)
    {
        correction = (gMonoMass[i] * gMultiplier * 10) - (gMonoMass_x100[i] * 10);
        if (correction >= 0)
        {
            gNodeCorrection[i] = correction + 0.5;
        }
        else
        {
            gNodeCorrection[i] = correction - 0.5;
        }
    }
    for (i = gAminoAcidNumber; i < MAX_GAPLIST; i++)
    {
        gNodeCorrection[i] = 0;
    }



/*      Create integer values for the monoisotopic masses of the six elements used by the program.*/

    for (i = 0; i < ELEMENT_NUMBER; i++)
    {
        gElementMass_x100[i] = (gElementMass[i] * gMultiplier) + 0.5;
    }

/*      Create integer values for gElementCorrection.*/

    for (i = 0; i < ELEMENT_NUMBER; i++)
    {
        correction = (gElementMass[i] * gMultiplier * 10) - (gElementMass_x100[i] * 10);
        if (correction >= 0)
        {
            gElementCorrection[i] = correction + 0.5;
        }
        else
        {
            gElementCorrection[i] = correction - 0.5;
        }
    }

/*      Create the appropriate integer values for the mass variables found in the .params file.*/

    gParam.peptideMW = gParam.peptideMW * gMultiplier;
    gParam.monoToAv = gParam.monoToAv * gMultiplier;
    gParam.peptideErr = gParam.peptideErr * gMultiplier;
    gParam.fragmentErr = gParam.fragmentErr * gMultiplier;
    gParam.qtofErr = gParam.qtofErr * gMultiplier;
    gParam.ionOffset = gParam.ionOffset * gMultiplier;
    gParam.cysMW = gParam.cysMW * gMultiplier;
    gParam.tagNMass = gParam.tagNMass * gMultiplier;
    gParam.tagCMass = gParam.tagCMass * gMultiplier;
    gParam.peakWidth = gParam.peakWidth * gMultiplier;
    gParam.modifiedNTerm = gParam.modifiedNTerm * gMultiplier;
    gParam.modifiedCTerm = gParam.modifiedCTerm * gMultiplier;

    gAvMonoTransition = AV_MONO_TRANSITION * gMultiplier;
    gWater = WATER * gMultiplier;
    gAmmonia = AMMONIA * gMultiplier;
    gCO = CO * gMultiplier;
    gAvResidueMass = AV_RESIDUE_MASS * gMultiplier;

/*      Convert the list of real data peaks to integer data peaks.*/

    currPtr = firstMassPtr;
    while (currPtr != NULL)
    {
        currPtr->mOverZ = (INT_4)(currPtr->mOverZ * gMultiplier + 0.5);
        currPtr = currPtr->next;
    }


    return;
}

/*
//--------------------------------------------------------------------------------
//  FreeAllSequenceScore()
//--------------------------------------------------------------------------------
    Free linked list of SequenceScore structs.
*/
void FreeSequenceScore(struct SequenceScore *currPtr)

{

    struct SequenceScore *freeMePtr;

    while (currPtr != NULL)
    {

        freeMePtr = currPtr;

        currPtr = currPtr->next;

        free(freeMePtr);

    }

    return;

}


/*
//--------------------------------------------------------------------------------
//  FreeMassList()
//--------------------------------------------------------------------------------
    Free linked list of MSData
*/
void FreeMassList(struct MSData *currPtr)

{

    struct MSData *freeMePtr;

    while (currPtr != NULL)
    {

        freeMePtr = currPtr;

        currPtr = currPtr->next;

        free((INT_4*)freeMePtr);                /****(INT_4*)****/

    }

    return;

}


/*
//--------------------------------------------------------------------------------
//  FreeSequence()
//--------------------------------------------------------------------------------
    Used for freeing memory in a linked list.  
*/
void FreeSequence(struct Sequence *currPtr)

{

    struct Sequence *freeMePtr;

    while (currPtr != NULL)
    {

        freeMePtr = currPtr;

        currPtr = currPtr->next;

        free(freeMePtr);

    }

    return;

}





/*
//--------------------------------------------------------------------------------
//  ReadDetailsFile()
//--------------------------------------------------------------------------------
    This function reads in the ion values from the Lutefisk.details editable ascii 
    file. The ion values are used to determine the values to assign to the nodes in
    the sequence graph that is to be developed later in the program. In addition to
    assigning values to each of the ion types to be used, there are three possible 
    ways of setting up the sequence graph.  The first is called "General", and is 
    an all-purpose fragmentation pattern where very little is assumed about the peptide.
    The second is called "Tryptic", and it assumes that the CID data is for a tryptic
    multiply charged precursor and that the CID was performed in a quadrupole 
    instrument under low energy CID conditions (this might also work for ion trap data
    - I don't know, I really just don't know).  The third type of fragmentation pattern
    is called "Arg+1", and it is for singly charged precursor ions that contain arginine.
    The type of fragmentation pattern is passed using the variable "fragmentPattern" 
    and is either 'G', 'T', or 'R', and depending on the value of "fragmentPattern" 
    either the first, second, or third column of Lutefisk.details is read to input 
    the various ion values.
*/
INT_4 ReadDetailsFile(void)
{
    FILE *fp;

    char stringBuffer[256];
    INT_4 totalIonVal = 0;
    INT_4 i = 0;
    INT_4 value;
    REAL_4 valueMultiplier;


    fp = fopen(gParam.detailsFilename, "r");

    if (fp == NULL)
    {
        printf("Cannot open Lutefisk details file.");
        exit(1);
    }

    while (!feof(fp))
    {
        if (my_fgets(stringBuffer, 256, fp) == NULL)
        {
            continue;
        }
        i += 1;

        if (gParam.fragmentPattern == 'G')
        {  /*   Read in ion values for the general fragmentation pattern.*/
            sscanf(stringBuffer, "%d %*d %*d", &value);
        }
        else if (gParam.fragmentPattern == 'T'  || gParam.fragmentPattern == 'Q')
        {  /*   Read in ion values for the triple quad tryptic fragmentation pattern.*/
            sscanf(stringBuffer, "%*d %d %*d", &value);
        }
        else if (gParam.fragmentPattern == 'L')
        {  /*   Read in ion values for the ion trap tryptic fragmenation pattern.*/
            sscanf(stringBuffer, "%*d %*d %d", &value);
        }


        if (i == 1)
        {
            gWeightedIonValues.b = (INT_4)value;
        }
        else if (i == 2)
        {
            gWeightedIonValues.a = (INT_4)value;
        }
        else if (i == 3)
        {
            gWeightedIonValues.c = (INT_4)value;
        }
        else if (i == 4)
        {
            gWeightedIonValues.d = (INT_4)value;
        }
        else if (i == 5)
        {
            gWeightedIonValues.b_minus17or18 = (INT_4)value;
        }
        else if (i == 6)
        {
            gWeightedIonValues.a_minus17or18 = (INT_4)value;
        }
        else if (i == 7)
        {
            gWeightedIonValues.y = (INT_4)value;
        }
        else if (i == 8)
        {
            gWeightedIonValues.y_minus2 = (INT_4)value;
        }
        else if (i == 9)
        {
            gWeightedIonValues.y_minus17or18 = (INT_4)value;
        }
        else if (i == 10)
        {
            gWeightedIonValues.x = (INT_4)value;
        }
        else if (i == 11)
        {
            gWeightedIonValues.z_plus1 = (INT_4)value;
        }
        else if (i == 12)
        {
            gWeightedIonValues.w = (INT_4)value;
        }
        else if (i == 13)
        {
            gWeightedIonValues.v = (INT_4)value;
        }
        else if (i == 14)
        {
            gWeightedIonValues.b_minusOH = (INT_4)value;
        }
        else if (i == 15)
        {
            gWeightedIonValues.b_minusOH_minus17 = (INT_4)value;
        }
    }

    fclose(fp);     

    totalIonVal =   gWeightedIonValues.b 
                    + gWeightedIonValues.a
                    + gWeightedIonValues.c
                    + gWeightedIonValues.d
                    + gWeightedIonValues.b_minus17or18
                    + gWeightedIonValues.a_minus17or18
                    + gWeightedIonValues.y
                    + gWeightedIonValues.y_minus2
                    + gWeightedIonValues.y_minus17or18
                    + gWeightedIonValues.x
                    + gWeightedIonValues.z_plus1
                    + gWeightedIonValues.w
                    +     gWeightedIonValues.v
                    + gWeightedIonValues.b_minusOH
                    + gWeightedIonValues.b_minusOH_minus17;

    if (totalIonVal > 30)    /*can't exceed value of a char, and 30 is well below 127*/
    {
        valueMultiplier = 30 / (REAL_4)totalIonVal;
        gWeightedIonValues.b                                    = (REAL_4)gWeightedIonValues.b * valueMultiplier + 0.5;
        gWeightedIonValues.a                                    = (REAL_4)gWeightedIonValues.a * valueMultiplier + 0.5;
        gWeightedIonValues.c                                    = (REAL_4)gWeightedIonValues.c * valueMultiplier + 0.5;
        gWeightedIonValues.d                                    = (REAL_4)gWeightedIonValues.d * valueMultiplier + 0.5;
        gWeightedIonValues.b_minus17or18                = (REAL_4)gWeightedIonValues.b_minus17or18 * valueMultiplier + 0.5;
        gWeightedIonValues.a_minus17or18                = (REAL_4)gWeightedIonValues.a_minus17or18 * valueMultiplier + 0.5;
        gWeightedIonValues.y                                    = (REAL_4)gWeightedIonValues.y * valueMultiplier + 0.5;
        gWeightedIonValues.y_minus2                     = (REAL_4)gWeightedIonValues.y_minus2 * valueMultiplier + 0.5;
        gWeightedIonValues.y_minus17or18                = (REAL_4)gWeightedIonValues.y_minus17or18 * valueMultiplier + 0.5;
        gWeightedIonValues.x                                    = (REAL_4)gWeightedIonValues.x * valueMultiplier + 0.5;
        gWeightedIonValues.z_plus1                              = (REAL_4)gWeightedIonValues.z_plus1 * valueMultiplier + 0.5;
        gWeightedIonValues.w                                    = (REAL_4)gWeightedIonValues.w * valueMultiplier + 0.5;
        gWeightedIonValues.v                                    = (REAL_4)gWeightedIonValues.v * valueMultiplier + 0.5;
        gWeightedIonValues.b_minusOH                    = (REAL_4)gWeightedIonValues.b_minusOH * valueMultiplier + 0.5;
        gWeightedIonValues.b_minusOH_minus17    = (REAL_4)gWeightedIonValues.b_minusOH_minus17 * valueMultiplier + 0.5;
        totalIonVal = 30;
    }

    return totalIonVal;

}

/*
//--------------------------------------------------------------------------------
//  ReadEdmanFile()
//--------------------------------------------------------------------------------
    This function reads the data from the file Lutefisk.edman into the INT_4 array 
    called gEdmanData[MAX_PEPTIDE_LENGTH][gAminoAcidNumber].  The array gEdmanData 
    contains the nominal masses of the amino acids listed in each cycle, and 
    gMaxCycleNum contains the number of cycles listed in lutefisk.edman.
*/
void ReadEdmanFile()
{

    FILE *fp;

    char stringBuffer[256], test;
    char edmanChar[MAX_PEPTIDE_LENGTH][AMINO_ACID_NUMBER];
    INT_4 i, j, k;

/* Initialize some variables.*/
    test = TRUE;
    gMaxCycleNum = 0;
    for (i = 0; i < MAX_PEPTIDE_LENGTH; i++)
    {
        for (j = 0; j < gAminoAcidNumber; j++)
        {
            gEdmanData[i][j] = 0;
            edmanChar[i][j] = 0;
        }
    }

/*Open the file.*/
    fp = fopen(gParam.edmanFilename, "r");

    if (fp == NULL)
    {
        printf("Could not open the Edman file '%s'. Quitting.", gParam.edmanFilename);
        exit(1);
    }

/*Read the information into a character array of single letter amino acid codes.*/
    while (!feof(fp) && test)
    {
        test = FALSE;
        if (my_fgets(stringBuffer, 256, fp) == NULL)
        {
            continue;
        }
        j = 0;
        while (stringBuffer[j] >= 65 && stringBuffer[j] <= 121 && j <= 19)
        {
            test = TRUE;
            if (stringBuffer[j] >= 97)
            {
                stringBuffer[j] = stringBuffer[j] - 32;
            }
            edmanChar[gMaxCycleNum][j] = stringBuffer[j];
            j+=1;
        }
        edmanChar[gMaxCycleNum][j] = 0;
        gMaxCycleNum += 1;
    }

    if (test == FALSE)       /*If it got out of the while loop by having test = FALSE, then the value of
                                             gMaxCycleNum will be one too many.*/
    {
        gMaxCycleNum = gMaxCycleNum - 1;
    }

    fclose(fp);

/*      Convert single letter code characters to nominal mass values.*/
    for (i = 0; i < gMaxCycleNum; i++)
    {
        j = 0;
        if (edmanChar[i][j] == 'X')
        {
            k = 0;
            for (j = 0; j < gAminoAcidNumber; j++)
            {
                if (gGapList[j] != 0)
                {
                    gEdmanData[i][k] = gGapList[j];
                    k++;
                }
            }
        }
        else
        {
            while (edmanChar[i][j] != 0)
            {
                for (k = 0; k < gAminoAcidNumber; k++)
                {
                    if (gSingAA[k] == edmanChar[i][j])
                    {
                        if (gGapList[k] == 0)
                        {
                            if (gSingAA[k] == 'Q' || gSingAA[k] == 'I')      /*If gGapList is 0, because its
                                                                                                                     Ile or Gln, then I need to make
                                                                                                                     special provision for these two
                                                                                                                     amino acids.  In the end, these
                                                                                                                     will come out looking like they
                                                                                                                     are Lys and Leu, so maybe I'll
                                                                                                                     need to fix that minor problem
                                                                                                                     later on.*/
                            {
                                if (gSingAA[k] == 'Q')
                                {
                                    gEdmanData[i][j] = gGapList[K];
                                }
                                if (gSingAA[k] == 'I')
                                {
                                    gEdmanData[i][j] = gGapList[L];
                                }
                            }
                            else
                            {
                                printf("An absent amino acid was found in the Edman data file.");
                                exit(1);
                            }
                        }
                        else
                        {
                            gEdmanData[i][j] = gGapList[k];
                        }
                    }
                }
                j++;
            }
        }
    }

    return;

}

/*
//--------------------------------------------------------------------------------
//  ReadParamsFile()
//--------------------------------------------------------------------------------
    This function reads parameters from the Lutefisk.params editable text file. 
    Fields set with command-line options override the values found in the 
    Lutefisk.params file. The values are stored in various fields of a global struct
    called gParam.  The meanings of the fields are as follows:

    fMonitor    = TRUE or FALSE and turns simulated teletype interface on/off (default = on)
    fVerbose    = TRUE or FALSE and spits out extra infor if on (default = off)
    paramFile   = string (default = lutefisk.params)
    outputFile  = string (default = lutefisk.out)
    detailsFilename = string (default = lutefisk.details)
    cidFilename = name of the asci file containing the CID data.
    peptideMW   = molecular weight of the peptide (Da).
    chargeState = charge on the precursor ion.
    fragmentErr = fragment ion mass tolerance (Da). 
    ionOffset   = any error in mass measurement that is consistent throughout the spectrum.
                  This value is added to the observed m/z values in cidFilename.
    fragmentPattern = G, T, L, or D for general, tryptic triple quad, tryptic LCQ 
                      fragmentations or default, where the program determines if its T or L
                      automatically.
    cysMW = molecular weight of cysteine, including any possible modification.
    proteolysis = T, K, E or D for tryptic, Lys-C, V8 digestion, or Asp-N.  This is used in 
                  setting up the graph.
    centroidOrProfile = C, P, or D depending on if the CID data is centroid or profile data.
                        D is the default value where the program determines this automatically.
    monoToAv = the mass (Da) above which the observed masses are assumed to be average mass
               and below which the observed masses are assumed to be monoisotopic.
               The cutoff is gradual, and occurs over a range extending 400 Da below monoToAv.
    ionsPerWindow = the number of ions in a SPECTRAL_WINDOW_WIDTH Da wide window, for 
                    multiply charged ions the window size becomes correspondingly narrower.
    aaPresent = a character string of amino acid single letter codes w/o spaces that represent
                the amino acids known to be present in the peptide.
    aaAbsent  = a character string of amino acid single letter codes w/o spaces that represent
                the amino acids known to be absent in the peptide.
    modifiedNTerm = N, A, C, or P for none, acetylated, carbamylated, or pyroglutamylated 
                    N-terminus.
    modifiedCterm = N or A for none or amidated C-terminal modification.
    tagNMass  = mass (Da) of the amino acid residues N-terminal to the sequence tag.
    tagSequence = a character string of amino acid single letter codes w/o spaces that represent
                  the amino acid sequence of the sequence tag.
    tagCMass  = mass (Da) of the amino acid residues C-terminal to the sequence tag.
    finalSeqNum = the number of completed sequences that will be saved for final scoring.
    topSeqNum = the number of subsequences allowed..
    extThresh = a fraction applied to the top scoring sequence extension.  All other extensions
                must have scores exceeding this value.
    maxExtNum = the maximum number of extensions allowed per subsequence.
    maxGapNum = the maximum number of two amino acid gaps in the subsequence.
    peakWidth = the width of the ions at 10% of the height.  A default of zero will cause
                the program to automatically determine the peak width if the data is in
                profile mode.  If the data is centroided, then it defaults to a value of 2.
    ionThreshold = the ion threshold times the average intensity in the spectrum is the
                   theshold below which signals are discarded.  The m/z's above the precursor
                   use a threshold that is one-half of ionThreshold.
    autoTag = Y or N.  Yes will initiate the automatic sequence tag finder.
    peptideErr = the peptide molecular weight tolerance (in Da).
    edmanDataFile = The filename where Edman data is located.
    ionsPerResidue = once the CID data is input, there is an upper limit on the total number of
                     ions that will be used in the analysis.  This upper limit is based on the
                     ionsPerResidue x the number of "average residues" in the peptide.  This
                     "average residue" number is determined from the AV_RESIDUE_MASS value
                     found in Lutefisk.h and peptideMW.
    CIDfileType = F or T, depending on if the file is ASCII generated by the Finnigan TSQ List
                  program, or a simple tab-delineated list.
*/
void  ReadParamsFile(void)
{
    FILE *   fp;
    char     stringBuffer[256], Leu = FALSE, Ile = FALSE, Gln = FALSE, Lys = FALSE;
    INT_4    i, j, length;
    char *   setting;
    char *   value;


    gParam.edmanPresent = FALSE;    /*set to true if an edman file is present*/

    fp = fopen(gParam.paramFile, "r");
    if (fp == NULL)
    {
        /* fopen() returns NULL if it can't find the file. */
        printf("Cannot open the parameter file '%s'\n", gParam.paramFile);
        goto problem;
    }

    while (!feof(fp))        /*While not at the end of file.*/
    {
        if (my_fgets(stringBuffer, 256, fp) == NULL)
        {
            /* fgets returns NULL at the end of a file.
               Drop out of the current iteration, so that the  
               while condition will terminate the loop.*/
            continue;       
        }

        setting = strtok(stringBuffer, ":");

        /* Skip comment lines */    
        if (setting[0] == '/') continue;

        value = strtok(NULL, " \t");
        if (NULL == value) continue;

        /* Skip settings w/o values. */
        if (value[0] == '|') continue;

/*        printf("Token: '%s' Value: '%s'\n", setting, value);    
*/
        if (!strcmp(setting, "CID Filename")) /*----------------------------------*/
        {
            if (!strlen(gParam.cidFilename))
            {
                strcpy(gParam.cidFilename, value);
            }
            if (gParam.fVerbose) printf("CID file name = %s\n", gParam.cidFilename);
        }
        else if (!strcmp(setting, "CID Quality"))
        {
            gParam.quality = toupper(value[0]);
            if (gParam.quality == 'Y')
            {
                gParam.quality = TRUE;
            }
            else
            {
                gParam.quality = FALSE;
            }
            if (gParam.fVerbose) printf("Quality = %d\n", gParam.quality);
        }
        else if (!strcmp(setting, "Peptide MW"))  /*------------------------------*/
        {
            gParam.peptideMW = atof(value);
            if (gParam.fVerbose) printf("peptide MW = %f\n", gParam.peptideMW);
        }
        else if (!strcmp(setting, "Charge-state"))  /*----------------------------*/
        {
            gParam.chargeState = atoi(value);
            if (gParam.fVerbose) printf("charge state = %d\n", gParam.chargeState);
        }
        else if (!strcmp(setting, "MaxEnt3"))  /*---------------------------------*/
        {
            gParam.maxent3 = toupper(value[0]);
            if (gParam.maxent3 == 'Y')
            {
                gParam.maxent3 = TRUE;
            }
            else
            {
                gParam.maxent3 = FALSE;
            }
            if (gParam.fVerbose) printf("MaxEnt3 = %d\n", gParam.maxent3);
        }
        else if (!strcmp(setting, "Peptide Error (u)"))  /*-----------------------*/
        {
            gParam.peptideErr = atof(value);
            if (gParam.peptideErr < 0 || gParam.peptideErr > 5)
            {
                printf("The peptide error should be positive and less than or equal to 5.\n");
                goto problem;
            }
            if (gParam.fVerbose) printf("peptide error = %f\n", gParam.peptideErr);
        }
        else if (!strcmp(setting, "Fragment Error (u)")) /*-----------------------*/
        {
            gParam.fragmentErr = atof(value);
            if (gParam.fragmentErr < 0 || gParam.fragmentErr > 5)
            {
                printf("The fragment ion error should be positive and less than or equal to 5.\n");
                goto problem;
            }
            if (gParam.fVerbose) printf("Fragment ion error = %f\n", gParam.fragmentErr);
        }
        else if (!strcmp(setting, "Final Fragment Err (u)"))  /*------------------*/
        {
            gParam.qtofErr = atof(value);
            if (gParam.qtofErr < 0 || gParam.qtofErr > 0.4)
            {
                printf("The qtof final fragment error should be less than 0.4.\n");
                goto problem;
            }
            if (gParam.fVerbose) printf("Qtof fragment error = %f\n", gParam.qtofErr);
        }
        else if(!strcmp(setting, "Score threshold"))
        {
        	gParam.outputThreshold = atof(value);
        	if(gParam.outputThreshold <= 0 || gParam.outputThreshold >= 1)
        	{
        		printf("The output score threshold should be between 0 and 1\n");
        		goto problem;
        	}
        	if(gParam.fVerbose) printf("Final output score threshold = %f\n", gParam.outputThreshold);
        }
        else if (!strcmp(setting, "Max. Final Sequences"))   /*------------------*/
        {
            gParam.finalSeqNum = atoi(value);
            if (gParam.finalSeqNum < 0)
            {
                printf("The number of completed sequences must be positive.\n");
                goto problem;
            }
            if (gParam.fVerbose) printf("Number of seqs to store = %d\n", gParam.finalSeqNum);
        }
        else if (!strcmp(setting, "Mass Scrambles for Statistics"))
        {
            gParam.wrongSeqNum = atoi(value);
            if (gParam.wrongSeqNum < 0)
            {
                printf("The number of mass scrambles is zero or higher.\n");
                goto problem;
            }
        }
        else if (!strcmp(setting, "Number of sequences"))
        {
        	gParam.outputSeqNum = atoi(value);
        	if(gParam.outputSeqNum <= 0 || gParam.outputSeqNum > 50)
        	{
        		printf("The number of output sequences must be between 0 and 50.\n");
        		goto problem;
        	}
        	if(gParam.fVerbose) printf("Number of output sequences to display in final report = %d\n", gParam.outputSeqNum);
        }
        else if(!strcmp(setting, "Shoe Size (US)"))
        {
        	gParam.shoeSize = atoi(value);
        	if(gParam.shoeSize > 15) printf("Your feet are enormous.\n");
        	if(gParam.shoeSize < 5)  printf("Your shoes are probably too tight.\n");
        }
        else if (!strcmp(setting, "Max. Subsequences")) /*----------------------*/
        {
            gParam.topSeqNum = atoi(value);
            if (gParam.topSeqNum < 0)
            {
                printf("The number of subsequences must be positive.\n");
                goto problem;
            }
            if (gParam.fVerbose) printf("Number of seqs to display = %d\n", gParam.topSeqNum);
        }
        else if (!strcmp(setting, "CID File Type"))  /*------------------------*/
        {
            gParam.CIDfileType = toupper(value[0]);

            /* If the input file ends w/ .dat then it's Native and if .dta it's a dta file 
               (unless it's been set as (Q) a Micromass .dta file which is a screwy variant). */
            length = strlen(gParam.cidFilename);
            if ((length > 4) && (!strncmp(gParam.cidFilename + length - 4, ".dta", 4))
                && gParam.CIDfileType != 'Q')
            {
                gParam.CIDfileType = 'D';
            }
            else if ((length > 4) && (!strncmp(gParam.cidFilename + length - 4, ".dat", 4)))
            {
                gParam.CIDfileType = 'N';
            }

            if (gParam.CIDfileType != 'F'       /* ICIS text file */
                && gParam.CIDfileType != 'T'    /* tab text file */
                && gParam.CIDfileType != 'L'    /* LCQ text file */
                && gParam.CIDfileType != 'N'    /* Finnigan '.dat' file */
                && gParam.CIDfileType != 'D'    /* Finnigan '.dta' file */
                && gParam.CIDfileType != 'Q')   /* Micromass pkl pseudo '.dta' format */
            {
                printf("Unregnized CID file type '%c'\n", gParam.CIDfileType);
                goto problem;
            }

            if (gParam.fVerbose) printf("CID file type = %c\n", gParam.CIDfileType);
        }
        else if (!strcmp(setting, "Profile/Centroid"))  /*--------------------*/
        {
            gParam.centroidOrProfile = toupper(value[0]);

            if (gParam.centroidOrProfile != 'C' 
                && gParam.centroidOrProfile != 'P'
                && gParam.centroidOrProfile != 'A')
            {
                printf("Illegal value '%c' for the data centroid or profile type.\n", gParam.centroidOrProfile);
                goto problem;
            }

            if (gParam.fVerbose) printf("Centroid/Profile = %c\n", gParam.centroidOrProfile);
        }
        else if (!strcmp(setting, "Peak Width (u)"))  /*----------------------*/
        {
            gParam.peakWidth = atof(value);

            gParam.peakWidth = gParam.peakWidth / 2;        
            /* Previously this was hard-coded as HALF_WINDOW,
               but it seemed like a good idea to make this a
               parameter in Lutefisk.gParam.  However, the concept
               of "peakwidth" was more intuitive than "half-window",
               hence the change.  To compensate, I divide the peak-
               width by two to get the half-window.  Thus, I can
               simply replace HALF_WINDOW in the old code with 
               peakWidth in the newer code without needing
               to worry about the differences.*/

            if (gParam.peakWidth < 0 || gParam.peakWidth > 5)
            {
                printf("Peakwidth of zero invokes the autopeak find.\n");
                printf("Otherwise choose a postive value less than or equal to 5\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Peak width = %f\n", gParam.peakWidth);
        }
        else if (!strcmp(setting, "Ion Threshold"))  /*-----------------------*/
        {
            gParam.ionThreshold = atof(value);

            if (gParam.ionThreshold < 0 || gParam.ionThreshold >= 10)
            {
                printf("The ion threshold should be positive and less than or equal to 10.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Ion threshold = %f\n", gParam.ionThreshold);
        }
        else if (!strcmp(setting, "Mass Offset (u)"))  /*---------------------*/
        {
            gParam.ionOffset = atof(value);

            if (gParam.ionOffset < -2 || gParam.ionOffset > 2)
            {
                printf("I think you're ion offset is kind of weird.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Ion offset = %f\n", gParam.ionOffset);
        }
        else if (!strcmp(setting, "Ions Per Window"))  /*---------------------*/
        {
            gParam.ionsPerWindow = atof(value);

            if (gParam.ionsPerWindow < 0)
            {
                printf("The ions per window should be positive.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Ions per window = %.1f\n", gParam.ionsPerWindow);
        }
        else if (!strcmp(setting, "Ions Per Residue"))  /*--------------------*/
        {
            gParam.ionsPerResidue = atof(value);

            if (gParam.ionsPerResidue < 0 || gParam.ionsPerResidue > 20)
            {
                printf("The ions per residue should be positive.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Ions per residue = %.1f\n", gParam.ionsPerResidue);
        }
        else if (!strcmp(setting, "Transition Mass (u)"))  /*-----------------*/
        {
            gParam.monoToAv = atoi(value);

            if (gParam.monoToAv < 0)
            {
                printf("The switch mass should be positive.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Switch mass = %d\n", gParam.monoToAv);
        }
        else if (!strcmp(setting, "Fragmentation Pattern"))   /*---------------*/
        {
            gParam.fragmentPattern = toupper(value[0]);

            if (gParam.fragmentPattern    != 'L'    /* Tryptic ion trap data */
                && gParam.fragmentPattern != 'G'        /* Currently G = general pattern
                                                           and is not really supported or recommended*/
                && gParam.fragmentPattern != 'T'        /* Tryptic triple quad data */
                && gParam.fragmentPattern != 'D'        /* Default that signals the program
                                                           to decide on its own whether the data
                                                           is from a triple quad or an ion trap*/
                && gParam.fragmentPattern != 'Q')       /* QTOF data */
            {
                printf("Lutefisk.gParam:  triple quad, qtof or ion trap tryptic fragmentation pattern?\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Fragmentation pattern = %c\n", gParam.fragmentPattern);
        }
        else if (!strcmp(setting, "Max. Gaps"))   /*---------------------------*/
        {
            gParam.maxGapNum = atoi(value);

            if (gParam.maxGapNum < -1 || gParam.maxGapNum > 5)
            {
                printf("The number of gaps per sequence should be less than 5.\n");
                printf("A value of -1 signals an automatic gap determination based on peptide mass.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Max gaps = %d\n", gParam.maxGapNum);
        }
        else if (!strcmp(setting, "Extension Threshold"))   /*----------------*/
        {
            gParam.extThresh = atof(value);

            if (gParam.extThresh < 0 || gParam.extThresh >= 1)
            {
                printf("The extension threshold is a value from 0 to 1.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Extension threshold = %f\n", gParam.extThresh);
        }
        else if (!strcmp(setting, "Max. Extensions"))   /*--------------------*/
        {
            gParam.maxExtNum = atoi(value);

            if (gParam.maxExtNum < 0 || gParam.maxExtNum > 10)
            {
                printf("The number of extensions is a positive number less than 10.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Max extensions = %d\n", gParam.maxExtNum);
        }
        else if (!strcmp(setting, "Cysteine Mass"))   /*----------------------*/
        {
            gParam.cysMW = atof(value);

            if (gParam.cysMW < 0)
            {
                printf("Lutefisk.gParam: The cysteine molecular weight must be positive.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Cysteine MW = %f\n", gParam.cysMW);
        }
        else if (!strcmp(setting, "Proteolysis"))   /*------------------------*/
        {
            gParam.proteolysis = toupper(value[0]);

            if (gParam.proteolysis    != 'T' 
                && gParam.proteolysis != 'K' 
                && gParam.proteolysis != 'E'
                && gParam.proteolysis != 'D'
                && gParam.proteolysis != 'N')
            {
                printf("Lutefisk.gParam: Type of proteolysis must be specified as:\n"
                       "tryptic (T), Lys-c (K), V8 (E), Asp-N (D), or none (N).\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Proteolysis = %c\n", gParam.proteolysis);
        }
        else if (!strcmp(setting, "Modified N-terminus"))    /*---------------*/
        {
            gParam.modifiedNTerm = atof(value);
            
            if (gParam.modifiedNTerm < 0 
                || (gParam.modifiedNTerm > gParam.peptideMW 
                    && gParam.peptideMW > 0))
            {
                printf("The N-terminal mass is unreasonable.  It is the R in R-NH-\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Modified N-terminus = %f\n", gParam.modifiedNTerm);
        }
        else if (!strcmp(setting, "Modified C-terminus"))    /*---------------*/
        {
            gParam.modifiedCTerm = atof(value);

            if (gParam.modifiedCTerm < 0 
                || (gParam.modifiedCTerm > gParam.peptideMW 
                    && gParam.peptideMW > 0))
            {
                printf("The C-terminal mass is unreasonable.  It is the R in -CO-R\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Modified C-terminus = %f\n", gParam.modifiedCTerm);
        }
        else if (!strcmp(setting, "Auto Tag"))    /*--------------------------*/
        {
            gParam.autoTag = toupper(value[0]);

            if (gParam.autoTag == 'Y')
            {
                gParam.autoTag = TRUE;
            }
            else
            {
                gParam.autoTag = FALSE;
            }

            if (gParam.fVerbose) printf("Auto tag? = %d\n", gParam.autoTag);
        }
        else if (!strcmp(setting, "Tag Low Mass y Ion"))    /*----------------*/
        {
            gParam.tagCMass = atof(value);

            if (gParam.tagCMass < 0 
                || (gParam.tagCMass > gParam.peptideMW 
                    && gParam.peptideMW > 0))
            {
                printf("The sequence tag N-terminal mass is unreasonable.\n");
                goto problem;
            }

            if (gParam.fVerbose) printf("Low mass y ion = %f\n", gParam.tagCMass);
        }
        else if (!strcmp(setting, "Sequence Tag"))   /*-----------------------*/
        {
            strcpy(gParam.tagSequence, value);

            for (i = 0; i < strlen(gParam.tagSequence); i++)
            {
                gParam.tagSequence[i] = toupper(gParam.tagSequence[i]);

                /* If Q or I are entered as part of the tag, then change them to the isobaric K and L.*/
                if (gParam.tagSequence[i] == 'Q'
                    && gParam.fragmentErr > gMonoMass[K] - gMonoMass[Q])
                {
                    gParam.tagSequence[i] = 'K';
                }

                if (gParam.tagSequence[i] == 'I')
                {
                    gParam.tagSequence[i] = 'L';
                }

            }

            if (gParam.fVerbose) printf("Sequence tag = %s\n", gParam.tagSequence);
        }
        else if (!strcmp(setting, "Tag High Mass y Ion"))   /*----------------*/
        {
            gParam.tagNMass = atof(value);

            if (gParam.tagNMass < gParam.tagCMass 
                || (gParam.peptideMW > 0 
                    && gParam.tagNMass > gParam.peptideMW + gElementMass[HYDROGEN] + gParam.fragmentErr))
            {
                if (gParam.tagSequence[0] != '*')
                {
                    printf("The sequence tag C-terminal mass is unreasonable.\n");
                    goto problem;
                }
            }

            if (gParam.fVerbose) printf("Low high y ion = %f\n", gParam.tagNMass);
        }
        else if (!strcmp(setting, "Present Amino Acids"))   /*----------------*/
        {
            strcpy(gParam.aaPresent, value);

            /* Force to uppercase */
            for (i = 0; i < strlen(gParam.aaPresent); i++)
            {
                gParam.aaPresent[i] = toupper(gParam.aaPresent[i]);
            }

            if (strchr(gParam.aaPresent, 'I')
                && !strchr(gParam.aaPresent, 'L'))
            {
                strcat(gParam.aaPresent, "L");
            }
            else if (strchr(gParam.aaPresent, 'L')
                     && !strchr(gParam.aaPresent, 'I'))
            {
                /*		strcat(gParam.aaPresent, "I");*/
            }

            if (gParam.fragmentErr > 0.4)
            {
                if (strchr(gParam.aaPresent, 'K')
                    && !strchr(gParam.aaPresent, 'Q'))
                {
                    strcat(gParam.aaPresent, "Q");
                }
                else if (strchr(gParam.aaPresent, 'Q')
                         && !strchr(gParam.aaPresent, 'K'))
                {
                    strcat(gParam.aaPresent, "K");
                }
            }

            if (gParam.fVerbose) printf("Amino acids present = %s\n", gParam.aaPresent);
        }
        else if (!strcmp(setting, "Absent Amino Acids"))    /*----------------*/
        {
            strcpy(gParam.aaAbsent, value);

            /* Force to uppercase */
            for (i = 0; i < strlen(gParam.aaAbsent); i++)
            {
                gParam.aaAbsent[i] = toupper(gParam.aaAbsent[i]);
            }

            if (strchr(gParam.aaAbsent, 'I')
                && !strchr(gParam.aaAbsent, 'L'))
            {
                strcat(gParam.aaAbsent, "L");
            }
            else if (strchr(gParam.aaAbsent, 'L')
                     && !strchr(gParam.aaAbsent, 'I'))
            {
                strcat(gParam.aaAbsent, "I");
            }

            if (gParam.fragmentErr > 0.4)
            {
                if (strchr(gParam.aaAbsent, 'K')
                    && !strchr(gParam.aaAbsent, 'Q'))
                {
                    strcat(gParam.aaAbsent, "Q");
                }
                else if (strchr(gParam.aaAbsent, 'Q')
                         && !strchr(gParam.aaAbsent, 'K'))
                {
                    strcat(gParam.aaAbsent, "K");
                }
            }

            if (gParam.fVerbose) printf("Amino acids absent = %s\n", gParam.aaAbsent);
        }
        else if (!strcmp(setting, "Edman Data File"))
        {
/*
*           The data is read into an array gEdmanData[cycle number][amino acids in the cycle], and
*           it contains the nominal mass values of the amino acids found in the specified file 
*           (rather than a character listing of the amino acid single letter code).
*/
            gParam.edmanPresent = TRUE;
            strcpy(gParam.edmanFilename, value);
            if (gParam.fVerbose) printf("Edman present = %d\n", gParam.edmanPresent);
            if (gParam.fVerbose) printf("Edman file name = %s\n", gParam.edmanFilename);
        }
        else if (!strcmp(setting, "DB Sequence File"))
        {
            strcpy(gParam.databaseSequences, value);
            if (gParam.fVerbose) printf("Database sequence file name = %s\n", gParam.databaseSequences);
        }
        else
        {
            printf("Unrecognized token '%s' in %s.\n", setting, gParam.paramFile);
            goto problem;
        }


    }         

    fclose(fp);


/*  Do multi-field interaction checks. */


/*      
*   Make sure that the amino acids that are absent (aaAbsent) do not match with the ones that
*   are present (aaPresent) or in the sequence tag (tagSequence).
*/
    if (gParam.aaAbsent[0] != '*')
    {
        if (gParam.aaPresent[0] != '*')
        {
            i = 0;
            while (gParam.aaAbsent[i] != 0)
            {
                j = 0;
                while (gParam.aaPresent[j] != 0)
                {
                    if (gParam.aaPresent[j] == gParam.aaAbsent[i])
                    {
                        printf("Amino acid '%c' has been listed as both present and absent.", gParam.aaPresent[j]);
                        goto problem;
                    }
                    j++;
                }
                i++;
            }
        }

        if (gParam.tagSequence[0] != '*')
        {
            i = 0;
            while (gParam.aaAbsent[i] != 0)
            {
                j = 0;
                while (gParam.tagSequence[j] != 0)
                {
                    if (gParam.tagSequence[j] == gParam.aaAbsent[i])
                    {
                        printf("An amino acid has been listed as a sequence tag and absent.");
                        goto problem;
                    }
                    j++;
                }
                i++;
            }
        }
    }

/*  Turn autotag off if a single tag is entered */
    if ((gParam.tagSequence[0] != '*' && gParam.tagCMass != 0 && gParam.tagNMass != 0) 
        && gParam.autoTag == TRUE)
    {
        gParam.autoTag = FALSE; 
    }


    if ((gParam.CIDfileType == 'D' || gParam.CIDfileType == 'X')
        && gParam.centroidOrProfile != 'C')
    {
        printf("Forcing .dta file to be read as centroid data.\n");
        gParam.centroidOrProfile = 'C'; /*force it to read centroid for dta files*/
    }

/*  To keep the qtof final scoring in sinc with the fragment error*/
    if (gParam.fragmentErr > MULTIPLIER_SWITCH / 10)
    {
        gParam.qtofErr = 0;
    }
    if (gParam.qtofErr >= gParam.fragmentErr)
    {
        gParam.qtofErr = 0;
    }
    if (gParam.fragmentPattern != 'Q')
    {
        gParam.qtofErr = 0;
    }

/*  Make sure maxent3 processing only allowed for qtof data*/
    if (gParam.maxent3)
    {
        if (gParam.fragmentPattern != 'Q')
        {
            gParam.maxent3 = FALSE;
        }
    }

/*	Make sure that if the peak width is zero, that no auto peakfinding is done for centroided
    data, and instead give reasonable values for different instruments*/
    if (gParam.centroidOrProfile == 'C' && gParam.peakWidth == 0)
    {
        if (gParam.fragmentPattern == 'T')
        {
            gParam.peakWidth = 1.5;
        }
        else if (gParam.fragmentPattern == 'Q')
        {
            gParam.peakWidth = 0.375;
        }
        else
        {
            gParam.peakWidth = 0.5;
        }
    }

/*	Check that the number of wrong sequences is reasonable, and then guess as to what was meant.*/
    if (gParam.wrongSeqNum < 3)
    {
        gParam.wrongSeqNum = 0; /*neg nums no good, and 1 or 2 is not statistically significant*/
    }
    /*For odd numbers, round up*/
    gParam.wrongSeqNum = (REAL_4)gParam.wrongSeqNum / 2 + 0.5;
    gParam.wrongSeqNum = gParam.wrongSeqNum * 2;


    return;

    problem:

    printf("\nQuitting.\n");
    exit(1);        
}

/*
//--------------------------------------------------------------------------------
//  ReadResidueFile()
//--------------------------------------------------------------------------------
    This file reads the amino acid residue masses.  Up to 25 residues are possible, including the
    twenty common ones.  Although Q/K have the same nominal masses, they can be differentiated
    if the error tolerances are less than about 0.04 u.  Although I and L are isomeric, I use the space
    for L to hold Leu or Ile and use the I position later in the program to store the values for
    oxidized Met, which can be differentiated from Phe if the error is less than 0.03 or so. 
    Modifications to Cys can be made via an entry in Lutefisk.params.  This leaves 5 positions --
    J, O, U, X, and Z for additional modified amino acids.
*/
void     ReadResidueFile(void)
{
    FILE *fp;

    char stringBuffer[256];
    char singleAA;
    REAL_4 monoisotopic, average;
    INT_4 nominal;
    INT_4 i = 0;
    INT_4 monoToNom;


    fp = fopen(gParam.residuesFilename, "r");

    if (fp == NULL)
    {
        printf("Cannot open Lutefisk residues file.");
        exit(1);
    }

    gAminoAcidNumber = 0;

    while (!feof(fp))
    {
        if (my_fgets(stringBuffer, 256, fp) == NULL)
        {
            continue;
        }

        sscanf(stringBuffer, "%c %f %f %d", &singleAA, &monoisotopic, &average, &nominal);

        if (monoisotopic != 0 && average != 0 && nominal != 0)
        {
            gSingAA[gAminoAcidNumber] = singleAA;
            gMonoMass[gAminoAcidNumber] = monoisotopic;
            gAvMass[gAminoAcidNumber] = average;
            gNomMass[gAminoAcidNumber] = nominal;
            gAminoAcidNumber++;
        }
    }

    /*Check for mistakes*/
    for (i = 0; i < gAminoAcidNumber; i++)
    {
        if (gMonoMass[i] < gAvMass[i] - 1 || gMonoMass[i] > gAvMass[i] + 1)
        {
            printf("*************************************************\n");
            printf("You may have a typo in Lutefisk.residues for %c.\n", gSingAA[i]);
            printf("*************************************************\n");
        }
        monoToNom = gMonoMass[i];
        if (monoToNom != gNomMass[i])
        {
            printf("*************************************************\n");
            printf("You may have a typo in Lutefisk.residues for %c.\n", gSingAA[i]);
            printf("*************************************************\n");
        }
    }

    return;
}


/*
//--------------------------------------------------------------------------------
//  SetupGapList()
//--------------------------------------------------------------------------------
    This function uses the values of cysMW and gNomMass to obtain a list of one and 
    two amino acid extensions.  Positions 0-19 contain the standard single amino acid
    extensions, where the values correspond to the nominal residue mass of the amino 
    acid. Subsequent index values contain the two amino acid extensions.  If any amino
    acids are listed in aaAbsent, then these are not used to generate either the two or
    one amino acid extensions. The number of values in the array gGapList is variable
    and depends on the number of amino acids that are absent - amino acids that are 
    absent are not used to make two amino acid extensions.  The number of extensions 
    in the array gGapList is gGapListIndex.
*/
void SetupGapList()
{
    INT_4 i, j, k, sum;
    char delAmAcid;
    INT_4 absentFlag;
    INT_4 duplicateFlag;
    INT_4 massToAdd;
    INT_4 lysGlnDiff = (gMonoMass_x100[K] - gMonoMass_x100[Q]) * 1.5;       /*mass diff between Q and K plus 50%*/
    struct MSData *currPtr = NULL;
    struct MSData *nextPtr = NULL;

    /*     char singleAAPresent[AMINO_ACID_NUMBER];
           REAL_4 massDiff; 
    */

    gGapListIndex = -1;

    for (i = 0; i < gAminoAcidNumber; i++)   /*Copy the single aa extension masses.*/
    {
        absentFlag = FALSE;
        if (gParam.aaAbsent[0] != '*') /* Check to see if the AA is on the absent list. */
        {
            /* (We won't add it to the gap list if it is.)   */
            delAmAcid = gParam.aaAbsent[0]; 
            j = 0;
            while (delAmAcid != 0 && (delAmAcid >= 'A' && delAmAcid <= 'Y'))
            {
                if (gSingAA[i] == delAmAcid)
                {
                    absentFlag = TRUE;
                    break;
                }
                j++;
                delAmAcid = gParam.aaAbsent[j];
            }
        }

        if (absentFlag || i == I || (i == Q && gParam.fragmentErr >= lysGlnDiff)) /* Ile and Gln, which are represented by Leu and Lys.*/
        {
            massToAdd = 0;
        }
        else if (i == C && gParam.cysMW != 0)
        {
            massToAdd = (gParam.cysMW) + 0.5;       /*Change the mass for cysteine (in case its alkylated).*/
        }
        else
        {
            massToAdd = gMonoMass_x100[i];
        }

        gGapList[i] = massToAdd;
    }       

    gGapListIndex = gAminoAcidNumber - 1;
    for (i = 0; i < gAminoAcidNumber; i++)   /*Fill in the masses of the 2 AA extensions.*/
    {
        for (j = i; j < gAminoAcidNumber; j++)
        {
            if (gGapList[i] == 0 || gGapList[j] == 0) continue;
            sum = gGapList[i] + gGapList[j];
            /*sum = ((gMonoMass[i] + gMonoMass[j]) * gMultiplier) + 0.5;*/

            duplicateFlag = FALSE;
            for (k = 0; k <= gGapListIndex; k++)
            {
                if (gGapList[k] > sum - gParam.fragmentErr 
                    && gGapList[k] < sum + gParam.fragmentErr)
                {

                    /* We already have this mass so don't add it to the list. */
                    duplicateFlag = TRUE;
                    break;
                }
            }

            /*for(k = 0; k < gAminoAcidNumber; k++)
            {
                    if(gGapList[k] != 0)
                    {
                            for(m = gAminoAcidNumber; m <= gGapListIndex; m++)
                            {
                                    if(sum - gGapList[k] <= gGapList[m] + gParam.fragmentErr
                                            && sum - gGapList[k] >= gGapList[m] - gParam.fragmentErr)
                                    {
                                            duplicateFlag = true;
                                            break;
                                    }
                            }
                    }
            }*/

            if (!duplicateFlag)
            {
                gGapListIndex = gGapListIndex + 1;
                gGapList[gGapListIndex] = sum;
            }
        }
    }

    /*get rid of single amino acids that aren't there (save time later)*/
    /*    if(firstMassPtr == NULL || firstMassPtr->next == NULL)
        {
            printf("LutefiskMain: firstMassPtr = NULL");
            exit(1);
        }
        for(i = 0; i < gAminoAcidNumber; i++)
        {
            singleAAPresent[i] = 0;	
        }
        currPtr = firstMassPtr;
        while(currPtr->next != NULL)
        {
            nextPtr = currPtr->next;
            while(nextPtr != NULL)
            {
                massDiff = nextPtr->mOverZ - currPtr->mOverZ; 
                if(massDiff <= gMonoMass_x100[W] + gParam.fragmentErr
                    && massDiff >= gMonoMass_x100[G] - gParam.fragmentErr )	
                {
                    
                    for(i = 0; i < gAminoAcidNumber; i++)
                    {
                        if(gMonoMass_x100[i] != 0)
                        {
                            if(massDiff <= gMonoMass_x100[i] + gParam.fragmentErr &&
                                massDiff >= gMonoMass_x100[i] - gParam.fragmentErr)
                            {
                                singleAAPresent[i] = 1;
                                break;
                            }
                        }
                    }
                }
                nextPtr = nextPtr->next;
            }
            currPtr = currPtr->next;
        }
        for(i = 0; i < gAminoAcidNumber; i++)
        {
            gGapList[i] = gGapList[i] * singleAAPresent[i];
        }*/

    return;
}

/*
//--------------------------------------------------------------------------------
//  SetupSequenceTag()
//--------------------------------------------------------------------------------
    Convert from highY, lowY, tag to gParam.tagNMass, gParam.tagCMass, and
    gParam.tagSequence.  Initially I had the user enter the unsequenced N-terminal
    mass, unsequenced C-terminal mass, and the sequence tag; however, in order to 
    simplify the use of tags, I now mimic the PeptideSearch protocol where the
    low mass y ion is entered as a float, the high mass y ion is entered as a float,
    and the sequence tag is entered from low mass to high mass, which is C-terminal
    to N-terminal direction.  For practical reasons (I've never used a b ion sequence
    tag), I only consider y ion tags.
*/
void SetupSequenceTag()
{
    char tag[256];
    REAL_4 highY, lowY, tagMass, aToMFactor;
    INT_4 i, k, j;

    /*Save old values for recalculations*/
    highY = gParam.tagNMass;
    lowY = gParam.tagCMass;
    i = 0;
    while (gParam.tagSequence[i] != 0)
    {
        tag[i] = gParam.tagSequence[i];
        i++;
    }
    tag[i] = 0;

    if (tag[0] == '*' || highY == 0 || lowY == 0)
    {
        gParam.tagSequence[0] = '*';    /*Make sure all of the tag paramters are set 
                                                                        properly if there is no tag.*/
        gParam.tagSequence[1] = '0';
        gParam.tagCMass = 0;
        gParam.tagNMass = 0;
    }
    else
    {
        /*Make some quick checks to see if I should shit-can the process*/
        if (highY > gParam.peptideMW)
        {
            printf("Your high mass y ion of your sequence tag is greater than the peptide mass.");
            exit(1);
        }
        if (highY < 0)
        {
            printf("Curiously, you have chosen a negative mass value for the high mass y ion.");
            exit(1);
        }
        if (lowY > gParam.peptideMW)
        {
            printf("Your low mass y ion of your sequence tag is greater than the peptide mass.");
            exit(1);
        }
        if (lowY < 0)
        {
            printf("Curiously, you have chosen a negative mass value for the low mass y ion.");
            exit(1);
        }
        printf("Specified sequence tag:  [%8.3f] %s [%8.3f]\n", lowY, tag, highY);
        gParam.tagCMass = lowY - (2 * gElementMass[HYDROGEN]);
        gParam.tagNMass = gParam.peptideMW - (highY - (2 * gElementMass[HYDROGEN]));
        j = 0;
        k = 0;
        while (tag[k] != 0)
        {
            k++;
        }
        k--;
        while (tag[j] != 0)
        {
            gParam.tagSequence[j] = tag[k];
            k--;
            j++;
        }
    }

    if (gParam.tagSequence[0] != '*')        /*Make sure tag matches peptide MW.*/
    {
        tagMass = gParam.tagCMass + gParam.tagNMass;    
        j = 0;
        while (gParam.tagSequence[j] != 0)
        {
            for (k = 0; k < gAminoAcidNumber; k++)
            {
                if (gSingAA[k] == gParam.tagSequence[j])
                {
                    if (k == C)
                    {
                        tagMass = tagMass + gParam.cysMW;
                    }
                    else
                    {
                        tagMass = tagMass + gMonoMass[k];
                    }
                    break;
                }
            }
            j++;
        }
        if (tagMass > gParam.monoToAv)
        {
            tagMass = tagMass * MONO_TO_AV;
        }
        else if (tagMass > gParam.monoToAv - AV_MONO_TRANSITION)
        {
            aToMFactor = (tagMass - (gParam.monoToAv -  AV_MONO_TRANSITION)) / 
                         AV_MONO_TRANSITION;
            aToMFactor = (MONO_TO_AV - 1) * aToMFactor;
            aToMFactor = 1 + aToMFactor;
            tagMass = tagMass * aToMFactor;
        }
        if ((tagMass > gParam.peptideMW + gParam.peptideErr) ||
            (tagMass < gParam.peptideMW - gParam.peptideErr))
        {
            printf("Your peptide tag does not match the peptide molecular weight.\n");
            exit(1);
        }
    }
    return;
}

/*
//--------------------------------------------------------------------------------
//  SystemCheck()
//--------------------------------------------------------------------------------
*/
BOOLEAN SystemCheck(void)
{
    BOOLEAN testValue = TRUE;
    BOOLEAN big_endian;

#if defined __BIG_ENDIAN
    big_endian = TRUE;
#else
    big_endian = FALSE;
#endif
    {
        UINT_4 test[2] = {0x41424344, 0x0}; /* ASCII "ABCD" in big endian */
        /* printf ("%s\n", (char *)&test); */

        if (big_endian)
        {
            if (strcmp((char *)&test, "ABCD"))
            {
                printf("Program should be set to __LITTLE_ENDIAN in LutefiskDefinitions.h\n");
                testValue = FALSE;
            }
        }
        else if (!strcmp((char *)&test, "ABCD"))
        {
            printf("Program should be set to __BIG_ENDIAN in LutefiskDefinitions.h\n");
            testValue = FALSE;      
        }
    }

    if (sizeof(REAL_4) != 4)
    {
        printf("REAL_4  is %d bytes instead of 4.\n", sizeof(REAL_4));
        testValue = FALSE;
    }
    if (sizeof(REAL_8) != 8)
    {
        printf("REAL_8  is %d bytes instead of 8.\n", sizeof(REAL_8));
        testValue = FALSE;
    }
    if (sizeof(INT_2) != 2)
    {
        printf("INT_2  is %d bytes instead of 2.\n", sizeof(INT_2));
        testValue = FALSE;
    }
    if (sizeof(UINT_2) != 2)
    {
        printf("UINT_2  is %d bytes instead of 2.\n", sizeof(UINT_2));
        testValue = FALSE;
    }
    if (sizeof(INT_4) != 4)
    {
        printf("INT_4  is %d bytes instead of 4.\n", sizeof(INT_4));
        testValue = FALSE;
    }
    if (sizeof(UINT_4) != 4)
    {
        printf("UINT_4  is %d bytes instead of 4.\n", sizeof(UINT_4));
        testValue = FALSE;
    }
    if (sizeof(CHAR) != 1)
    {
        printf("CHAR  is %d bytes instead of 1.\n", sizeof(CHAR));
        testValue = FALSE;
    }
    if (sizeof(BOOLEAN) != 1)
    {
        printf("BOOLEAN  is %d bytes instead of 1.\n", sizeof(BOOLEAN));
        testValue = FALSE;
    }


    return(testValue);
}

/****************************PrintHeaderToFile***************************************************
*	This function prints header information to the output file.
*/

void PrintPartingGiftToFile()
{
    FILE *fp;

    PrintHeaderToFile();

    /* Open the file for appending.*/
    fp = fopen(gParam.outputFile, "a");
    if(fp == NULL)	/*fopen returns NULL if there's a problem.*/
    {
        printf("Cannot open %s for appending.\n", gParam.outputFile);
        exit(1);
    }

    printf("\n\nNo potential candidate sequences could be found.\n\n");
    fprintf(fp, "No potential candidate sequences could be found.\n\n");

    fclose(fp);

}

/****************************PrintHeaderToFile***************************************************
*	This function prints header information to the output file.
*/

void PrintHeaderToFile()
{
	INT_4 i;
	FILE *fp;
   	const 	time_t		theTime = (const time_t)time(NULL);
	

        /* Open a new file.*/
	fp = fopen(gParam.outputFile, "w");
	if(fp == NULL)	/*fopen returns NULL if there's a problem.*/
	{
		printf("Cannot open %s to write the output.\n", gParam.outputFile);
		exit(1);
	}

	fprintf(fp, versionString);
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
	 
        /* Pad the output with a few blank lines for later use.*/
	fprintf(fp, "\n\n\n");

	fclose(fp);
	return;
}


/******************************AdjustPeptideMW****************************************
*
*	Ion trap data usually has complementary pairs of ions (b and y ions for the same 
*	cleavage.  These pairs can be used to adjust the peptide MW, which is often not 
*	very accurately determined in the MS scan.
*/

void	AdjustPeptideMW(struct MSData *firstMassPtr)
{
	REAL_4	mass[MAX_ION_NUM], peptide[MAX_ION_NUM], testMass, maxPeptideMass;
	REAL_4	pairMass[MAX_ION_NUM], mass2[MAX_ION_NUM];
	REAL_4	lowMassIon[MAX_ION_NUM], water, ammonia, avePairMass;
	REAL_8	stDev;
	INT_4 	ionNum, i, j, peptideNum, nominalPeptide[MAX_ION_NUM], maxPeptideNum;
	INT_4 nominalCount, requiredPairs, pairNum;
	BOOLEAN test;
	struct MSData *currPtr;
	
	/*initialize*/
	water 			= gElementMass[OXYGEN] + 2 * gElementMass[HYDROGEN];
	ammonia 		= gElementMass[NITROGEN] + 3 * gElementMass[HYDROGEN];
	peptideNum 		= 0;
	ionNum 			= 0;
	maxPeptideNum 	= 0;
	maxPeptideMass	= 0;
	pairNum			= 0;
	stDev			= 0;
	avePairMass		= 0;
	
	for(i = 0; i < MAX_ION_NUM; i++)
	{
		peptide[i] 			= 0;
		nominalPeptide[i]	= 0;
		pairMass[i] 		= 0;
		mass[i] 			= 0;
		mass2[i]			= 0;
	}
	
	if(gParam.peptideMW < 750)
	{
		requiredPairs = 1;	/*need more than this number of pairs*/
	}
	else if(gParam.peptideMW < 1500)
	{
		requiredPairs = 2;
	}
	else if(gParam.peptideMW < 2250)
	{
		requiredPairs = 3;
	}
	else
	{
		requiredPairs = 4;
	}
	
/*    Fill in the mass array assuming singly-charged ions.*/
	currPtr = firstMassPtr;
	while(currPtr != NULL)
	{
		mass[ionNum] = currPtr->mOverZ - gElementMass[HYDROGEN];
		ionNum++;
		currPtr = currPtr->next;
	}
	
	/*    Fill in the mass array assuming doubly-charged ions.*/
	if(gParam.chargeState > 2)
	{
		for(i = 0; i < ionNum; i++)
		{
   			testMass = (mass[i] + gElementMass[HYDROGEN]) * 2 - 2 * gElementMass[HYDROGEN];
   			if(testMass < gParam.peptideMW - gMonoMass[G] + gParam.fragmentErr &&
   				testMass > 700)	/*doubly charged ions have to be in the right mass range*/
   			{
   		    	mass2[i] = testMass;
   		    }
   		    else
   		    {
   		    	mass2[i] = 0;	/*zero is a flag that it could not be a doubly-charged ion*/
   		    }
   		}
	}
	
	/*	Find a suitable error*/
	/*assume all ions are singly-charged*/
	for(i = 0; i < ionNum - 1; i++)
	{
		for(j = i + 1; j < ionNum; j++)
		{
			testMass = mass[i] + mass[j];
			if(testMass <= gParam.peptideMW + gParam.peptideErr * 2 &&
				testMass >= gParam.peptideMW - gParam.peptideErr * 2)
			{
				pairMass[pairNum] = testMass;	/*collect the data*/
				pairNum++;
			}
		}
	}
	/*now assume that one of the pair is doubly-charged*/
	if(gParam.chargeState > 2)
	{
		for(i = 0; i < ionNum - 1; i++)
		{
			for(j = i + 1; j < ionNum; j++)
			{
				if(mass2[j] > mass[i])
				{
					testMass = mass[i] + mass2[j];
					if(testMass <= gParam.peptideMW + gParam.peptideErr * 2 &&
						testMass >= gParam.peptideMW - gParam.peptideErr * 2)
					{
						pairMass[pairNum] = testMass;	/*collect the data*/
						pairNum++;
					}
				}
			}
		}
	}
	
	
	if(pairNum < 3)
	{
		stDev = gParam.peptideErr;	/*not enough pairs of ions to determine standard deviation
									so it gets defined as the peptide error from the params file*/
	}
	else	/*enough data to take a stab at finding standard deviation*/
	{
		for(i = 0; i < pairNum; i++)
		{
			avePairMass += pairMass[i];
		}
		avePairMass = avePairMass / pairNum;
		
		for(i = 0; i < pairNum; i++)
		{
			stDev += ((pairMass[i] - avePairMass) * (pairMass[i] - avePairMass));
		}
		stDev = stDev / (pairNum - 1);
		stDev = sqrt(stDev);
	
	}
	
	/*reality checks*/
	if(stDev > 2 * gParam.peptideErr)
	{
		stDev = 2 * gParam.peptideErr;	/*don't let the error be too big*/
	}
	else if(stDev < 0.5 * gParam.peptideErr)
	{
		stDev = 0.5 * gParam.peptideErr;	/*or too small*/
	}
	

	/*find pairs of masses that are close to the peptide molecular weight*/
	/*first assume the ions are all singly-charged*/
	for(i = 0; i < ionNum - 1; i++)
	{
		for(j = i + 1; j < ionNum; j++)
		{
			testMass = mass[i] + mass[j];
			if(testMass <= gParam.peptideMW + stDev &&
				testMass >= gParam.peptideMW - stDev)
			{
				peptide[peptideNum] = testMass;
				lowMassIon[peptideNum] = mass[i];
				peptideNum++;
			}
		}
	}
	/*now assume that one of them is doubly-charged*/
	if(gParam.chargeState > 2)
	{
		for(i = 0; i < ionNum - 1; i++)
		{
			for(j = i + 1; j < ionNum; j++)
			{
				if(mass2[j] > mass[i])
				{
					testMass = mass[i] + mass2[j];
					if(testMass <= gParam.peptideMW + stDev &&
						testMass >= gParam.peptideMW - stDev)
					{
						peptide[peptideNum] = testMass;
						lowMassIon[peptideNum] = mass[i];
						peptideNum++;
					}
				}
			}
		}
	}
	
	/*Find the correct nominal masses of the peptides that were identified*/
	for(i = 0; i < peptideNum; i++)
	{
		nominalPeptide[i] = peptide[i] - (peptide[i] * 0.00050275) + 0.5;
	}
	
	/*Wipe out nominalPeptides that are derived from ions that are too close in mass */
	for(i = 0; i < peptideNum; i++)
	{
		testMass = nominalPeptide[i];
		if(testMass != 0)
		{
			for(j = 0; j < peptideNum; j++)
			{
				if(nominalPeptide[j] == testMass && i != j)
				{
					
					if(fabs(lowMassIon[i] - lowMassIon[j]) <= gElementMass[HYDROGEN] + gParam.fragmentErr)
					{
						test = TRUE;
					}
					else
					{
						test = FALSE;
					}
					
				/*	test = TRUE;
				//	if(fabs(lowMassIon[i] - lowMassIon[j]) >= ammonia - gParam.fragmentErr &&
				//		fabs(lowMassIon[i] - lowMassIon[j]) <= ammonia + gParam.fragmentErr)
				//	{
				//		test = FALSE;	
				//	}
				//	if(fabs(lowMassIon[i] - lowMassIon[j]) >= water - gParam.fragmentErr &&
				//		fabs(lowMassIon[i] - lowMassIon[j]) <= water + gParam.fragmentErr)
				//	{
				//		test = FALSE;	
				//	}
				//	if(fabs(lowMassIon[i] - lowMassIon[j]) >= gMonoMass[G] - gParam.fragmentErr)
				//	{
				//		test = FALSE;	
				//	}*/
					if(test)	/*if the two ions differ by a small mass thats not water or ammonia*/
					{
						nominalPeptide[j] = 0;
					}
				}
			}
		}
	}
	
	/*Count the numbers of each nominal peptide mass*/
	for(i = 0; i < peptideNum; i++)
	{
		testMass = nominalPeptide[i];
		if(testMass != 0)
		{
			nominalCount = 0;
			for(j = 0; j < peptideNum; j++)
			{
				if(nominalPeptide[j] == testMass)
				{
					nominalCount++;
				}
			}
			if(nominalCount > maxPeptideNum)
			{
				maxPeptideNum = nominalCount;
				maxPeptideMass = testMass;
			}
		}
	}
	
	/*decide what to do with this information*/
	if(maxPeptideNum > requiredPairs)	/*need at least 3 pairs of ions to change the peptide MW*/
	{
		maxPeptideMass = maxPeptideMass + (maxPeptideMass * 0.00050275);	/*add estimated mass defect*/
		testMass = fabs(maxPeptideMass - gParam.peptideMW);
		if(testMass <= gParam.peptideErr * 1.5)
		{
			maxPeptideNum = maxPeptideNum - requiredPairs;	/*maxPeptideNum is now the number in excess of requiredPairs*/
			testMass = maxPeptideMass * maxPeptideNum + gParam.peptideMW;	
			testMass = testMass / (maxPeptideNum + 1);	/*calculate the average of the theoretical and obsv'd MW*/
			printf("Peptide MW was adjusted from %f to ", gParam.peptideMW);
			gParam.peptideMW = testMass;
			printf("%f using %ld ion pairs\n", gParam.peptideMW, maxPeptideNum + requiredPairs);
		}
	}
	

	return;
}













