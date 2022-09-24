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

/* ANSI Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Lutefisk Headers */
#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"
#include "ListRoutines.h"

/* Function prototypes */
tMSDataList *ReadFinniganFile(char *filename);
tMSDataList *ReadFinniganFile(char *filename) {
        printf("That would be nice, wouldn't it.\nNative file reading not implimented.\n");
        return NULL;
        exit(1);
}
/*    Declare pointers to a struct that is global only within LutefiskGetCID.c.    */
static struct MSData *gGroupingPtr;
static struct MSData *gLastDataPtr;

/*Definitions for this file*/
#define MIN_NUM_IONS 5  /*Minimum number of ions after processing in GetCID*/
#define MAX_ION_MASS 3000       /*Ions greater than this are deemed too high to not be a mistake*/
#define MIN_HIGHMASS_INT_RATIO 0.05     /*Ratio of high mass intensity over total intensity*/
#define HIGH_MASS_RATIO 0.95    /*Ions are counted until this % of high mass ion intensity is reached*/
#define LCQ_INT_SUM_CUTOFF 500  /*Cutoff for good intensity total for LCQ data*/
#define QTOF_INT_SUM_CUTOFF     140  /*Cutoff for good intensity total for Qtof data*/
#define MAX_HIGH_MASS 100       /*Max number of ions greater than precursor*/
#define MAX_MASS 4000   /*Peptides above this mass are tossed out.*/
#define MIN_MASS 800    /*Peptides below this mass are tossed out.*/
#define LOW_MASS_ION_NUM 19     /*Number of peptide-related low mass ions*/


/*
//--------------------------------------------------------------------------------
//  GetCidData()
//--------------------------------------------------------------------------------
    GetCidData uses an ASCII input file named cidFilename, where cidFilename is produced 
    using the 'Print...' command in the Finnigan program LIST.  When cidFilename is 
    placed in the same folder as Lutefisk, it finds weighted average masses within a 
    resolution window defined by the variable 'peakWidth'.  'peakWidth' is half of the 
    peak width at its base, and for high sensitivity (low resolution) applications peakWidth 
    is typically 4-6 Da, but is best set to 3.0 for unit resolution (this helps get rid 
    of the C13 peaks).  GetCidData returns the value of the MSData struct called firstAvMassPtr,
    which is the final list of ion m/z and intensity values used for the rest of the program.

    Modified 03.13.00 JAT - Split out the file reading code into ReadCIDFile().
*/

struct MSData *GetCidData(void)
{
    tMSDataList     *MSDataList = NULL;
    tMSDataList     *peakList   = NULL;
    struct MSData     *firstAvMassPtr = NULL;
    INT_4 i, finalIonCount;
    REAL_4    excessIonRatio;
    REAL_4 generalQuality = 0;
    REAL_4 lowMassQuality = 0;
    REAL_4 highMassQuality = 0;
    REAL_4 quality = 0;


    if(gParam.fMonitor)
    {
        printf("Processing CID datafile '%s'\n", gParam.cidFilename);
    }

    MSDataList = ReadCIDFile(gParam.cidFilename);

    if (MSDataList->numObjects == 0)
    {
        printf("There doesn't seem to be any data in the firstDataPtr linked list.\n");
        exit(1);
    }    
    
    TrimList(MSDataList);
 
    
/*
*    Find the global multiplier, which is the value that determines the number of decimal points
*    for which the data will be considered significant.  For example, if two decimal points are 
*    to be used, then the multiplier will be 100 (ie, 100 x the monoisotopic mass will be the 
*    node values used in the graph.
*/

    FindTheMultiplier();


/*
*    If the 'Autodetect' value for the centroidOrProfile variable in the .params file was,
*    selected, the program tries to figure out what type of file is being used here.
*/
    if(gParam.centroidOrProfile == 'A')
    {
        CentroidOrProfile(MSDataList);
    }

    
/*
*    If the default value for the fragmentation pattern was selected, then the program
*    decides if the data is TSQ or LCQ data.  The differentiation is not fool-proof, but
*    should work most of the time.  The idea is that LCQ msms data starts at a mass that is
*    1/3 of the precursor ion m/z.
*/
    if(gParam.fragmentPattern == 'D')
    {
        GuessAtTheFragmentPattern();
    }

    
/*
*    Here's where profile data is smoothed once with a 5 point smoothing routine 
*    (using Finnigan's coefficients). 
*/

    if(gParam.centroidOrProfile == 'P')
    {
        SmoothCID(MSDataList);
    }

       
/*    Here is where the auto-peakWidth is determined.  The lowest mass ions (less than 500,
    and no more than ten and most intense) are chosen from the raw data in firstDataPtr.  
    The peak tops in firstDataPtr are found, and the peak width at 50% are calculated.  
    The peakWidth is 2x this 50% peakWidth, since
    peakWidth assumes a 10% valley.  Peaks used must be >10% of the most intense ion below 500.
    An average of these peakWidths is determined, and outliers are discarded and the average
    is recalculated.  This recalculated values is the peakWidth in the automated mode.*/

    if(gParam.peakWidth == 0)
    {
        if(gParam.centroidOrProfile == 'P')
        {
            gParam.peakWidth = GetPeakWidth(MSDataList);
        }
        else
        {
            if(gParam.fragmentPattern == 'L')
            {
                gParam.peakWidth = 0.5;    /*for lcq, assume width of 1 da*/
            }
            else 
            {
                gParam.peakWidth = 1;    /*otherwise, assume width of 2 da*/
            }
        }
        printf("The peak width is %5.2f\n", gParam.peakWidth * 2);
    }


/* 
*    Find the signal threshold, where the noise is the average of all ions and the 
*    signal/noise threshold is defined by ionThreshold.
*/

    gParam.intThreshold = FindThreshold(MSDataList);
    
/*    Alternatively, you could use the median rather than an average.*/

    /*threshold = FindMedian(firstDataPtr);    TOO SLOW */

/* 
*    Identify groups of ions, where a group of ions are above the "threshold" and 
*    consecutive ions are less than a peakWidth apart in m/z value (this is of particular
*    concern when importing centroided data rather than profile data).  Groups of ions 
*    are placed in a separate linked list which is passed to the function 'IonSorter'.   
*    The global struct pointer gGroupingPtr is used to remember the position in the linked 
*    list of the complete CID data where IonGrouper is to start finding a new group.  If 
*    IonGrouper is no longer able to find ion groups (ie, end of list), it passes
*   a NULL value.  This signals IonSorter to also pass a NULL value, which signals 
*    the while loop to terminate ('test' becomes FALSE).  Otherwise, IonSorter returns a pointer
*    to a linked list of MSData structs, which contains a short list of ions that have been
*    weight averaged.  This list of ions is appended onto the linked list that is pointed to by
*    firstAvMassPtr.
*/

        
    peakList = IonCondenser(MSDataList);

    if(peakList->numObjects < 5)
    {
        printf("Too few data points remaining after ion condensation!");
        return NULL;
    }
    
/*    
*    After generating all of the mass averaged ion values, some ions may still be closer than 
*    the value of 'peakWidth'.  This happens because some ions are just outside of the peakWidth
*    tolerance, but after weighted averaging to determine their new masses, they get too close.
*    For pairs of ions that are too close together, those with the lowest intensity are removed.
*/

    /* XXXXXXXXXXX JAT - Is this still necessary? */
/*    firstAvMassPtr = ZeroTheIons(firstAvMassPtr);    
*/


/*    Eliminate ions below the threshold.*/
/*    currPtr = firstAvMassPtr->next;
    previousPtr = firstAvMassPtr;
    while(currPtr != NULL)
    {
        if(currPtr->mOverZ > precursor + gParam.fragmentErr)
            break;
        if(currPtr->intensity == threshold)
        {
            if((currPtr->mOverZ < 147.5) ||
                (currPtr->mOverZ < 175.5 && currPtr->mOverZ > 174.5) ||
                (currPtr->mOverZ < 159.5 && currPtr->mOverZ > 158.5))
            {
                currPtr = currPtr->next;
                previousPtr = previousPtr->next;
            }
            else
            {
                previousPtr->next = currPtr->next;
                free(currPtr);
                currPtr = previousPtr->next;
            }
        }
        else
        {
            currPtr = currPtr->next;
            previousPtr = previousPtr->next;
        }
    }
*/    


    
/*    
*    Next the summed intensity of the linked list of ions is checked to see if it is less
*    than 2 billion.  If the sum exceeds this then each ion is attenuated ten fold.  This
*    is to make sure that I don't try to use a number too big for a INT_4 to hold later
*    on.
*/

    CheckTheIntensity(peakList);

    
    
/*    Next add the ion offset to the m/z values in the linked list of ions.*/

    AddTheIonOffset(peakList);
    
/*    Remove low mass ions that are not due to amino acids.*/

    LowMassIonRemoval(peakList);
    
/*    For Qtof data, convert any intense doubly charged ions to singly-charged ones.*/

/*    
*    Next the program checks to see if ions that are 1 Da apart are due to isotopes.  The
*    ions that seem to be due to isotopes are removed.  This is not done for data processed
*        using maxent3, since the data would already have been de-isotoped.
*/

        if(!gParam.maxent3)
        {
                RemoveIsotopes(peakList);
        }


/*
*    Remove the precursor ions.
*/
    RemovePrecursors(peakList);
    
    
/*
*    Find the ions that might be y ions that can be connected via single amino acids to the
*    y1 ions of 147 and 175 (if they are present).  A value of 1 is placed in normIntensity
*    if the ion is a golden boy (from CA).
*/
    if(gParam.proteolysis == 'T' && (gParam.fragmentPattern == 'Q' || gParam.fragmentPattern == 'T'))
    {
        FindTheGoldenBoys(peakList);
    }
    
    if(gParam.fragmentPattern =='L' && peakList->numObjects < 100)
    {
        FindBYGoldenBoys(peakList);     /*spare the b/y pairs if there are not many ions at this point*/
    }
        
    
/*    
*    Next the program checks to see if there are too many ions clustered together.  It does
*    this by counting the number of ions within windows of width SPECTRAL_WINDOW_WIDTH Da and 
*    making sure that only a certain number of ions (ionsPerWindow) are present within any 
*    given window.  If there are too many ions, it throws out those with the lowest intensity.
*/

    WindowFilter(peakList);
    
/*
*    Find high mass ions that could not be either b or y ions.
*/

    EliminateBadHighMassIons(peakList);
    
/*
*    Verify that the selected ions have a signal to noise ratio greater than SIGNAL_NOISE that
*    is #defined in LutefiskDefinitions.  The dta files from qtof data are also checked for s/n problems.
*/
    if(gParam.centroidOrProfile == 'P' || gParam.CIDfileType == 'X' || gParam.CIDfileType == 'D')
    {
        CheckSignalToNoise(peakList, MSDataList);
    }
    
/*
*    Check to see if ions can be connected to other ions with very high mass accuracy 
*    (compared to the calculated amino acid masses).  For LCQ data I've decided (from 50 
*    measurements) that the average error should be better than 0.15.  For QTOF data, this
*    should be much higher.  For triple quad data, I don't care to speculate.
*/
    if(gParam.fragmentPattern == 'L' || gParam.fragmentPattern == 'Q')
    {
        CheckConnections(peakList);
    }
    
/*
*       For LCQ data, find b/y pairs and flag them as being golden boys that are difficult to delete on 
*       the basis of ion intensity.
*/
    
    if(gParam.fragmentPattern =='L')
    {
        FindBYGoldenBoys(peakList);
    }
  
/*    
*    Next I count the ions, and if there are more ions per residue than stipulated by the value
*    of "ionsPerResiude" (from Lutefisk.params), then the lowest intensity ions are removed.  
*    For example if the peptide MW is 1200 and there are supposed to be 5 ions per residue,
*    then only 50 of the most intense ions are kept.
*/
    
    finalIonCount = (gParam.peptideMW / AV_RESIDUE_MASS) + 0.5;    /*Num of average residues.*/
    finalIonCount = finalIonCount * gParam.ionsPerResidue;        /*Num of ions allowed.*/
    
    /*lowMassIons = 0;
    currPtr = &peakList->mass[0];
    ptrOfNoReturn = &peakList->mass[peakList->numObjects];
    while (currPtr < ptrOfNoReturn) 
    {
        if (currPtr->mOverZ < 146.5) 
        {
            lowMassIons++;
        }
        currPtr++;
    }
    finalIonCount += lowMassIons;*/
    
    if (peakList->numObjects  > finalIonCount)
    {
        WeedTheIons(peakList, finalIonCount, TRUE);
    }

/*
    The use of "golden boy" ions (ions that connect to 147 and 175 in tryptic peptides) can 
    produce too many ions.  The following procedure is done in order to reduce excessive numbers.
    The golden boys are lumped in with the regular ions here, ie, no special treatment.
*/
    excessIonRatio = (2000 - gParam.peptideMW) / 1000;
    excessIonRatio = (excessIonRatio * 0.25) + 1;
    if(excessIonRatio < 1)
        excessIonRatio = 1;

    finalIonCount = finalIonCount * excessIonRatio;    /*Allow 15% increase in number of ions.*/
    if (peakList->numObjects > finalIonCount)
    {
        WeedTheIons(peakList, finalIonCount, FALSE);
    }




/*    Print out the number of ions in the final linked list of averaged ions.*/

    if(gParam.fMonitor)
    {
        printf("Number of ions: %ld \n", peakList->numObjects);
    }

    if(peakList->numObjects < 5)
    {
        printf("Too few data points remaining after preprocessing!");
        return NULL;
    }
      
  



/*
*    For QTof data, use y1 ions for R and K and immonium ions to obtain a mass offset correction.
*/

    /*if(gParam.fragmentPattern == 'Q')
    {
        CalibrationCorrection(peakList);
    }*/

/*
*    For precursor charges of two or less, most of the fragment ions will be singly-charged.
*    Thus, one can compare theoretical mass defects from the observed defects, and make 
*    corrections.
*/

    if(gParam.chargeState <= 2 && gParam.fragmentPattern != 'Q')
    {
        DefectCorrection(peakList);
    }
    

/*
*    For LCQ data, normalize the intensity to the fourth most intense ions.  It seems that
*    the ion trap will produce one or two whopper ions, so it seems reasonable to avoid
*    having a large percentage of ion current in one or two ions.
*/

    if(gParam.fragmentPattern == 'L')
    {
        NormalizeIntensity(peakList);
    }

/*    Now get rid of the raw CID data.*/    
    DisposeList(MSDataList);
    
/*    Print the final list*/
    if(gParam.fVerbose)
    {
        DumpMSData(peakList);
    }


    /* KLUDGE:  Dump back into linked list - JAT */
/*    for (i = 0; i < MSDataList->numObjects; i++) 
    {    
        firstDataPtr = AddToCIDList(firstDataPtr, 
                                      LoadMSDataStruct(MSDataList->mass[i].mOverZ, 
                                                       MSDataList->mass[i].intensity));
    }        
*/    for (i = 0; i < peakList->numObjects; i++) 
    {    
        firstAvMassPtr = AddToCIDList(firstAvMassPtr, 
                                      LoadMSDataStruct(peakList->mass[i].mOverZ, 
                                                       peakList->mass[i].intensity));
    }        

/*      Check the CID data quality here, but only if there is a monitor for output.*/

        if(gParam.quality && gParam.fMonitor)
        {
                printf("\n");
                printf("Quality assessment:");
                
/*      Check the charge state, peptide mass, number of ions, total intensity of ions, and
        distribution of ions.*/
        
                generalQuality = GeneralEval(peakList);
        
/*      Check the low mass ions for qtof data*/

                if(gParam.fragmentPattern == 'Q')
                {
                        lowMassQuality = LowMassIonCheck(peakList);
                }
                else
                {
                        lowMassQuality = 1;     /*LCQ data doesn't contain the low mass end*/
                }

/*      Check for series of ions above the precursor that can be connected by amino acids*/

                highMassQuality = HighMassIonCheck(peakList);
        
/*      Combine the scores via multiplication; all values are between 0 and 1, so the combined
        quality value is also between 0 and 1*/
        
                quality = generalQuality * lowMassQuality * highMassQuality;
        
/*      Print some more info*/

                if(quality)
                {
                        printf("\nThis spectrum exceeds the minimal spectral quality parameters.\n");
                }
                else
                {
                        printf("\nThis spectrum stinks.\n");
                }
                /*printf("General appearance (scale of 0 to 1): %f\n", generalQuality);
                if(gParam.fragmentPattern != 'L')
                {
                        printf("Low mass ions (scale of 0 to 1): %f\n", lowMassQuality);
                }
                printf("High mass ions above the precursor (scale of 0 to 1): %f\n", highMassQuality);
                if(quality == 1)
                {
                        printf("Overall data quality rating (A to F scale): A\n\n");
                }
                else if (quality < 1 && quality > 0.8)
                {
                        printf("Overall data quality rating (A to F scale): B\n\n");
                }
                else if(quality <= 0.8 && quality >= 0.5)
                {
                        printf("Overall data quality rating (A to F scale): C\n\n");
                }
                else if(quality < 0.5 && quality > 0)
                {
                        printf("Overall data quality rating (A to F scale): D\n\n");
                }
                else
                {
                        printf("Overall data quality rating (A to F scale): F\n");
                        printf("Don't waste your time on this one.\n\n");
                        exit(0);
                }*/
                
        }

    DisposeList(peakList);
    
    return(firstAvMassPtr);
}

/*****************************GeneralEval******************************************************
*
*       Check that the charge state of the precursor is either +2 or +3, and that the mass is not
*       a weird value.  Also checks that this is not a virtually empty file (ie, more than a handful
*       of ions.  Checks that fragment ions are not weird values.  Makes sure that a minimum percentage
*       of the total ion abundance is above the precursor, and that of those above the precursor the
*       abundance is not isolated to just a few ions.  If all of these tests are passed, then a 
*       quality assignment is made based on the total intensity.
*
*/

REAL_4 GeneralEval(tMSDataList  *peakList)
{
        INT_4 i;
        INT_4 highMassIonCount = 0; 
        REAL_4 quality = 0;     
        REAL_4 precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN])
                                                / gParam.chargeState;
        REAL_4 totIntensity = 0;
        REAL_4 highMassInt = 0;
        REAL_4 highMassRatio = 0;
        REAL_4 highMassInt2 = 0;
        REAL_4 minPeakNum = gParam.peptideMW / AV_RESIDUE_MASS;
        tMSDataList *  intensityOrderedList = NULL;

/*      Do some basic checks on the charge state and peptide mass*/

        if(gParam.peptideMW < MIN_MASS)
        {
                return(0);      /*Too small*/
        }
        if(gParam.peptideMW > MAX_MASS)
        {
                return(0);      /*Too big*/
        }
        if(gParam.chargeState < 2)
        {
                return(0);      /*Not enough charge*/
        }
        if(gParam.chargeState > 3)
        {
                return(0);      /*Too much charge*/
        }

/*      Check the number of ions*/

        if(peakList->numObjects < minPeakNum)
        {
                return(0);      /*Not enough ions*/
        }
        if(peakList->numObjects < MIN_NUM_IONS)
        {
                return(0);
        }
        
/*      Check for weird ions (negative mass or very high)*/
        
        for (i = 0; i < peakList->numObjects; i++) 
    {    
        if(peakList->mass[i].mOverZ < -1)
        {
                return(0);
        }
        if(peakList->mass[i].mOverZ > MAX_ION_MASS)
        {
                return(0);
        }
    }
    
/*      Check for ions above the precursor*/
        
        for(i = 0; i < peakList->numObjects; i++)
        {
                totIntensity += peakList->mass[i].intensity;
                if(peakList->mass[i].mOverZ > precursor + gParam.fragmentErr * 3)
                {
                        highMassInt += peakList->mass[i].intensity;
                }
    }
    
    if(totIntensity == 0)
    {
        return(0);      /*No ion intensity, so quality of zero is returned*/
    }
    
/*      Calculate the percentage of total ion current above the precursor*/

    highMassRatio = highMassInt / totIntensity;
    
    if(highMassRatio < MIN_HIGHMASS_INT_RATIO)
    {
        return(0);      /*there's less than 5% of total intensity above the precursor*/
    }
    
/*      Now I count the number of ions that comprise 95% of the high mass intensity*/
    
    /*The next few lines are directly from bits of jat's code*/
    intensityOrderedList = (tMSDataList *) CopyList( peakList );
        if (!intensityOrderedList) 
    {
        printf("Ran out of memory in GeneralEval()!\n");
        exit(1);
    }
    /* Sort the intensityOrderedList in order of decreasing intensity */
    qsort(intensityOrderedList->mass,(size_t)intensityOrderedList->numObjects,
          (size_t)sizeof(tMSData),IntensityDescendSortFunc);
    
/*      Count the number of ions required to account for most of the high mass intensity.*/

    for(i = 0; i < intensityOrderedList->numObjects; i++)
    {
        if(intensityOrderedList->mass[i].mOverZ > precursor + 3 * gParam.fragmentErr)
        {
                highMassIonCount++;
                highMassInt2 += intensityOrderedList->mass[i].intensity;
        }
        if((highMassInt2 / totIntensity) > HIGH_MASS_RATIO * highMassRatio)
        {
                break;  /*For example if HIGH_MASS_RATIO = 0.9 and highMassRatio was 0.5 (ie,
                                half of the ion current in a spectrum was above the precursor), then
                                this bit of code counts the number of ions required to account for
                                0.45 of the total intensity (ie, 90% of the current above the precursor).
                                The idea is if only one or two ions account for most of the ion current,
                                then this is not as good as if the current were spread out among several
                                ions.  The latter situation would imply better sequence data, than just
                                having a few massive ions.*/
        }
    }
    
/*      Assign quality values based on the number of major high mass ions, and peptide mass.*/
    
    if(highMassIonCount == 1 && gParam.peptideMW > 850)
    {
        DisposeList(intensityOrderedList);      /*get rid of this ion list*/
        return(0);      /*Only one major ion above the precursor*/
    }
    
    DisposeList(intensityOrderedList);  /*get rid of this ion list*/
    
        return(1);      /*return a value of 1 if the spectrum passes this test*/
}

/*****************************LowMassIonCheck**************************************************
*
*       First, low mass ions are checked against an approved list of peptide-related ions.
*       If any are found, then the low mass quality is not zero, and might be as much as 1
*       depending on the number of low mass ions found.  If no peptide-related low mass ions
*       are found, then the number of weird ions are counted.  If any weird non-peptide ions
*       are located then the quality is zero.
*
*/

REAL_4 LowMassIonCheck(tMSDataList  *peakList)
{
        REAL_4 lowMassIons[LOW_MASS_ION_NUM] = {        60.0,   70.1,  72.1,    74.1,   84.0,   86.1,   
                                                                                                87.1,   88.0,   101.1,  102.1,  104.1,  110.1,  
                                                                                                112.1,  120.1,  129.1,  136.1,  159.1,  147.1,  
                                                                                                175.1
                              };
    REAL_4 quality = 1; /*default for no low mass ions, and no weird ions*/
        INT_4 i, j;
        INT_4 lowMassIonCount = 0;
        INT_4 weirdIons = 0;
        char test;
        
/*      First look for immonium and other peptide-related low mass ions*/

        for (i = 0; i < peakList->numObjects; i++)
        {
                if(peakList->mass[i].mOverZ > 175.1 + gParam.fragmentErr)
                {
                        break;
                }
                for(j = 0; j < LOW_MASS_ION_NUM; j++)
                {
                        if(peakList->mass[i].mOverZ <= lowMassIons[j] + gParam.fragmentErr &&
                                peakList->mass[i].mOverZ >= lowMassIons[j] - gParam.fragmentErr)
                        {
                                lowMassIonCount++;
                                break;
                        }
                }
        }
        
/*      If no peptide-related low mass ions found, then look for weird low mass ions.*/

        if(lowMassIonCount == 0)
        {
                for (i = 0; i < peakList->numObjects; i++)
                {
                        if(peakList->mass[i].mOverZ > 147.1 - gParam.fragmentErr)
                        {
                                break;
                        }
                        test = TRUE;    /*test is TRUE if its a weird peak*/
                        for(j = 0; j < 17; j++)
                        {
                                if(peakList->mass[i].mOverZ <= lowMassIons[j] + gParam.fragmentErr &&
                                        peakList->mass[i].mOverZ >= lowMassIons[j] - gParam.fragmentErr)
                                {
                                        test = FALSE;
                                        break;
                                }
                        }
                        if(test)
                        {
                                weirdIons++;
                        }
                }
        }
        
/*      Now assign a low mass ion quality.*/

        if(lowMassIonCount == 0)
        {
                if(weirdIons != 0)
                {
                        quality = 0;    /*Weird ions in the absence of immoniums is not good at all*/
                }
        }

        return(quality);
}


/*****************************HighMassIonCheck*************************************************
*
*
*
*/
REAL_4 HighMassIonCheck(tMSDataList  *peakList)
{
        INT_4 i, highCount, j, k, totalSinglyChargedConnectNum, totalDoublyChargedConnectNum;
        INT_4 nodeConnect[MAX_HIGH_MASS][AMINO_ACID_NUMBER];
        INT_4 connectNum[MAX_HIGH_MASS];
        INT_4 m, n, firstAA, secondAA, thirdAA, fourthAA;
        REAL_4 highMass[MAX_HIGH_MASS], highInt[MAX_HIGH_MASS];
        REAL_4 precursor;
        REAL_4 totalHighInt;
        REAL_4 averageInt, massDiff;
        REAL_4 threshold, quality;
        char runOfTwo, runOfThree, runOfFour;
        
/*      initialize*/
        
        averageInt = 0;
        totalHighInt = 0;
        quality = 0;
        precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN])
                                                / gParam.chargeState;
        highCount = 0;
        totalSinglyChargedConnectNum = 0;
        totalDoublyChargedConnectNum = 0;
        
        for(i = 0; i < MAX_HIGH_MASS; i++)
        {
                connectNum[i] = 0;
                for(j = 0; j < AMINO_ACID_NUMBER; j++)
                {
                        nodeConnect[i][j] = 0;
                }
        }
        for(i = 0; i < MAX_HIGH_MASS; i++)
        {
                highMass[i] = 0;
                highInt[i] = 0;
        }

        /*Identify high mass ions*/
        for (i = 0; i < peakList->numObjects; i++) 
        {       
                if(peakList->mass[i].mOverZ > precursor + 4 * gParam.fragmentErr)
                {
                        highMass[highCount] = peakList->mass[i].mOverZ;
                        highInt[highCount] = peakList->mass[i].intensity;
                        totalHighInt += peakList->mass[i].intensity;
                        highCount++;
                }
                if(highCount > MAX_HIGH_MASS)
                {
                        return(0);      /*Too many high mass ions (>100); I'll exceed the array sizes.
                                                Also, its unlikely to have this many high mass ions unless the
                                                file is screwed up somehow.*/
                }
        }
        
        if(highCount == 0)
        {
                return(0);      /*No high mass ions, so return a zero quality value*/
        }
        
        averageInt = totalHighInt / highCount;  /*average intensity*/
        threshold = averageInt / 5;                             /*threshold is 1/4 of the average; there is no
                                                                                        good reason to use 1/4 rather than, say 1/5.*/
        
        
        
/*      First assume fragment ions are all singly charged.  This first bit of nested for loops
        determines which nodes connect to each other.  The indexing for connectNum matches with
        highMass, and contains the number of connections that can be made to lower mass nodes.
        The array nodeConnect has two indexing, where one matches connectNum and highMass and the
        other index corresponds to the number of connections that can be made to lower mass nodes.
        The intensity of the higher mass node always has to be above the threshold, whereas the 
        lower mass node can be below threshold, but only if the mass difference corresponds to Pro
        or Gly.  The actual value stored by nodeConnect is the index value that it can connect with.*/

        for(i = 0; i < highCount; i++)
        {
                for(j = i; j < highCount; j++)
                {
                        if(highInt[j] > threshold)
                        {
                                massDiff = highMass[j] - highMass[i];
                                if(massDiff >= gMonoMass[G] - gParam.fragmentErr)
                                {
                                        for(k = 0; k < AMINO_ACID_NUMBER; k++)
                                        {
                                                if(massDiff <= gMonoMass[k] + gParam.fragmentErr &&
                                                        massDiff>= gMonoMass[k] - gParam.fragmentErr)
                                                {
                                                        if(highInt[i] > threshold || k == G || k == P)
                                                        {
                                                                nodeConnect[j][connectNum[j]] = i;
                                                                connectNum[j]++;
                                                                break;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
        
/*      Count the singly-charged connections; this number does not tell you the longest sequence 
        obtainable.  This is determined below....*/

        totalSinglyChargedConnectNum = 0;
        for(i = 0; i < highCount; i++)
        {
                if(connectNum[i] != 0)
                {
                        totalSinglyChargedConnectNum++;
                }
        }
        
/*      Can I make runs of two, three, or four amino acids (ie, can I link three, four, or five ions)?*/

        runOfTwo = FALSE;       /*two amino acids defined*/
        runOfThree = FALSE;     /*three amino acids defined*/
        runOfFour = FALSE;      /*four amino acids defined*/
        
        if(gParam.peptideMW > 1000)     /*smaller peptides might not give data with long series of nodes*/
        {
                for(i = highCount - 1; i >= 0; i--)     /*start at the high mass end and word down*/
                {
                        if(connectNum[i] != 0)  /*keep looking at the next node down if the current one
                                                                        lacks any connections to lower nodes*/
                        {
                                for(j = 0; j < connectNum[i]; j++)      /*loop through the all of the possible 
                                                                                                        connections*/
                                {
                                        firstAA = nodeConnect[i][j];    /*firstAA is the next node down*/
                                        if(connectNum[firstAA] == 0)
                                        {
                                                continue;       /*if the next node down does not connect to anything, then
                                                                        stop following this pathway*/
                                        }
                                        for(k = 0; k < connectNum[firstAA]; k++)
                                        {
                                                secondAA = nodeConnect[firstAA][k];     /*secondAA is the second node down*/
                                                runOfTwo = TRUE;        /*you can connect at least three ions, and this will
                                                                                        remain TRUE for the rest of this function*/
                                                if(connectNum[secondAA] == 0)
                                                {
                                                        continue;       /*if the second node down does not connect to anything, 
                                                                                then stop following this pathway*/
                                                }
                                                for(m = 0; m < connectNum[secondAA]; m++)       /*etc etc*/
                                                {
                                                        thirdAA = nodeConnect[secondAA][m];
                                                        runOfThree = TRUE;
                                                        if(connectNum[thirdAA] == 0)
                                                        {
                                                                continue;
                                                        }
                                                        for(n = 0; n < connectNum[thirdAA]; n++)
                                                        {
                                                                fourthAA = nodeConnect[thirdAA][n];
                                                                runOfFour = TRUE;
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }
                                                        
        
/*      Now assume fragment ions are doubly-charged, which is only possible if the precursor is 
        triply-charged, and the data was not derived from any maxent3 treatment (conversion of all
        ions to singly charged, and de-isotoped).*/
        
        if(gParam.chargeState == 3 && !gParam.maxent3)
        {
                /*initialize the node arrays to get rid of the singly charged info*/
                for(i = 0; i < MAX_HIGH_MASS; i++)
                {
                        connectNum[i] = 0;
                        for(j = 0; j < AMINO_ACID_NUMBER; j++)
                        {
                                nodeConnect[i][j] = 0;
                        }
                }
                
                for(i = 0; i < highCount; i++)
                {
                        for(j = i; j < highCount; j++)
                        {
                                if(highInt[j] > threshold)
                                {
                                        massDiff = (highMass[j] - highMass[i]) * 2;
                                        if(massDiff >= gMonoMass[G] - gParam.fragmentErr)
                                        {
                                                for(k = 0; k < AMINO_ACID_NUMBER; k++)
                                                {
                                                        if(massDiff <= gMonoMass[k] + gParam.fragmentErr &&
                                                                massDiff>= gMonoMass[k] - gParam.fragmentErr)
                                                        {
                                                                if(highInt[i] > threshold || k == G || k == P)
                                                                {
                                                                        nodeConnect[j][connectNum[j]] = i;
                                                                        connectNum[j]++;
                                                                        break;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
                /*Count the doubly-charged connections*/
                totalDoublyChargedConnectNum = 0;
                for(i = 0; i < highCount; i++)
                {
                        if(connectNum[i] != 0)
                        {
                                totalDoublyChargedConnectNum++;
                        }
                }
                
                /*Can I make runs of two, three, or four amino acids?*/
                
                if(gParam.peptideMW > 1000)
                {
                        for(i = highCount - 1; i >= 0; i--)
                        {
                                if(connectNum[i] != 0)
                                {
                                        for(j = 0; j < connectNum[i]; j++)
                                        {
                                                firstAA = nodeConnect[i][j];
                                                if(connectNum[firstAA] == 0)
                                                {
                                                        continue;
                                                }
                                                for(k = 0; k < connectNum[firstAA]; k++)
                                                {
                                                        secondAA = nodeConnect[firstAA][k];
                                                        runOfTwo = TRUE;
                                                        if(connectNum[secondAA] == 0)
                                                        {
                                                                continue;
                                                        }
                                                        for(m = 0; m < connectNum[secondAA]; m++)
                                                        {
                                                                thirdAA = nodeConnect[secondAA][m];
                                                                runOfThree = TRUE;
                                                                if(connectNum[thirdAA] == 0)
                                                                {
                                                                        continue;
                                                                }
                                                                for(n = 0; n < connectNum[thirdAA]; n++)
                                                                {
                                                                        fourthAA = nodeConnect[thirdAA][n];
                                                                        runOfFour = TRUE;
                                                                }
                                                        }
                                                }
                                        }
                                }
                        }
                }
        }

/*Print relevant info*/

        /*if(runOfFour)
        {
                printf("\n");
                printf("At least five ions above the precursor can be connected\n");
        }
        else if(runOfThree)
        {
                printf("\n");
                printf("At least four ions above the precursor can be connected\n");
        }
        else if(runOfTwo)
        {
                printf("\n");
                printf("At least three ions above the precursor can be connected\n");
        }
        else
        {
                if(totalSinglyChargedConnectNum == 0 && totalDoublyChargedConnectNum == 0)
                {
                        printf("\n");
                        printf("No ions above the precursor can be connected\n");
                }
                else
                {
                        printf("\n");
                        printf("Only pairs of ions above the precursor can be connected\n");
                }
        }*/
        
/*      First see if any connections can be made*/

        if(totalSinglyChargedConnectNum == 0 && totalDoublyChargedConnectNum == 0)
        {
                if(gParam.peptideMW < 950)
                {
                        quality = 1;    /*low mw peptides may not have many ions above the precursor*/
                }
                else
                {
                        quality = 0;    /*bigger peptides are therefore crap*/
                }
        }
        else    /*If any connections can be made, then start out w/ a high quality*/
        {
                quality = 1;
        }
        
/*      Quality is adjusted downward if insufficient sequence lengths obtained for certain mass ranges*/
        /*if(gParam.peptideMW > 1150)   /*Anything less than 1150 is not attenuated further*/
        /*{
                if(gParam.peptideMW < 1300)
                {
                        if(!runOfTwo)
                        {
                                quality *= 0.5;
                        }
                }
                else if(gParam.peptideMW < 1450)
                {
                        if(!runOfThree)
                        {
                                quality *= 0.5;
                        }
                }
                else
                {
                        if(!runOfFour)
                        {
                                quality *= 0.75;
                        }
                }
        }*/
        
        return(quality);
}

/*
//--------------------------------------------------------------------------------
//  ReadCIDFile()
//--------------------------------------------------------------------------------
    Modified 03.13.00 JAT - Split this function out from GetCidData().
*/
tMSDataList *ReadCIDFile(char *inFilename)
{

    FILE *     fp;
    INT_4      i;
    INT_4      j;
    tMSDataList * MSDataList   = NULL;
    tMSData    massToAdd;    
    REAL_4     massValue       = 0.0;
    REAL_4     oldMassValue    = 0.0;
    INT_4      ionIntensity    = 0;    
    INT_4      oldIonIntensity = 0;    
    REAL_4     intensityAsReal = 0.0;
    char *     stringBuffer    = NULL;
    char *     stringBuffer2   = NULL;
    BOOLEAN    firstNumberFlag = true;
    BOOLEAN    headerFlag      = true;

	msms.scanMassHigh = -1;
    stringBuffer = (char *)malloc(258);
    if (NULL == stringBuffer)
    {
        printf("Outa memory in GetCidData()!\n");
        goto problem;
    }
    
    stringBuffer2 = (char *)malloc(258);
    if (NULL == stringBuffer2)
    {
        printf("Outa memory in GetCidData()!\n");
        goto problem;
    }
    
    for (i = 0; i < 258; i++)
    {
        stringBuffer2[i] = 0;
    }

    if (gParam.CIDfileType == 'N')
    {
        /* Open the native data file and make a linked list of m/z and intensity values.*/
        if (gParam.fVerbose) printf("Reading the native CID file '%s'\n", inFilename);
        
        MSDataList = ReadFinniganFile(inFilename);
    }
    else 
    {
        /* Open the ASCII data file and make a linked list of m/z and intensity values.*/

        MSDataList = (tMSDataList *) CreateNewList( sizeof(tMSData), 500, 500 );
        if ( NULL == MSDataList )
        {
            printf("Outa memory in GetCIDFile()!\n");
            goto problem;
        }
    
        fp = fopen(inFilename,"r");
        if (fp == NULL)
        {
            printf("Cannot open the CID file '%s'.\n", inFilename);
            goto problem;
        }
    
        i=0;
        oldMassValue = 0;
        oldIonIntensity = 0;
    
        while (my_fgets(stringBuffer, 256, fp) != NULL)
        {    
            i+=1;
            
            /* Skip blank lines */
            if (!strcmp(stringBuffer, "\r")) continue;    /*PC files are screwy*/
            if (!strcmp(stringBuffer, "\n")) continue;    
        
            /* Deal with the headers first -------------------------------------- */
            if (headerFlag) 
            {
                if(gParam.CIDfileType == 'T')
                {
                        headerFlag = FALSE;     /*tab text has no header*/
                }
                else if (gParam.CIDfileType == 'F') /* Finnigan's ICIS text format */
                {
                    /* headerFlag is true until the data is being 
                    read w/ in the Finnigan file.*/
                    sscanf(stringBuffer, "%f %d", &massValue, &ionIntensity);
    
                    /* The number 1 should appear on the first 
                    line of data in a Finnigan ASCII file.*/
                    if(massValue == 1.0)
                    {
                        headerFlag = FALSE;
                    }
                    else if (massValue > 1.0 && massValue < 2000.0 && ionIntensity > 0)
                    {    /* Catch potential problem caused by forgetting
                           to change from 'F' to 'T' in the .params file. */
                        gParam.CIDfileType = 'T';
                    }
                    else continue; /* Still in the header, so go back to the start of the loop. */
                }        
                else if (gParam.CIDfileType == 'L') /* Finnigan's LCQ text format */
                {
                    sscanf(stringBuffer, "%[DataPeaks]", stringBuffer2);
                    if(!strcmp(stringBuffer2, "DataPeaks"))
                    {
                        headerFlag = FALSE;
                        continue;
                    }
                    else
                    {
                        continue;    /*read another line*/
                    }
                }
                else if (gParam.CIDfileType == 'D') /* Finnigan's '.dta' text format */
                {
                    /*First line is MH+ and charge: read if these
                    values are zero from the params file*/
                    if(gParam.peptideMW == 0 || gParam.chargeState == 0)
                    {
                        sscanf(stringBuffer, "%f %d", &gParam.peptideMW, &gParam.chargeState);
                        /*adjust for the fact that the value in the dta file is MH+*/
                        gParam.peptideMW -= gElementMass[HYDROGEN];
                        if (gParam.fVerbose) 
                        {
                           printf("  Precursor Mass:  %.3f\n", gParam.peptideMW);
                           printf("Precursor Charge:  %d\n", gParam.chargeState);
                        }
                    }
                    else
                    {    
                        sscanf(stringBuffer, "%*f %*d");
                    }
                    headerFlag = FALSE;
                    continue; /*the next line is the start of the data*/
                }
                else if (gParam.CIDfileType == 'Q') /* Micromass' pkl '.dta' text format */
                {
                    /*First line is the precursor m/z, followed by intensity (float) and charge: read if these
                    values are zero from the params file*/
                    if(gParam.peptideMW == 0 || gParam.chargeState == 0)
                    {
                        sscanf(stringBuffer, "%f %*f %d", &gParam.peptideMW, &gParam.chargeState);
                        /*adjust for the fact that the value in the pkl file is precurson mass*/
                        gParam.peptideMW = (gParam.peptideMW * gParam.chargeState) 
                                            - (gParam.chargeState * gElementMass[HYDROGEN]);
                    }
                    else
                    {    
                        sscanf(stringBuffer, "%*f %*f %*d");
                    }
                    headerFlag = FALSE;
                    continue; /*the next line is the start of the data*/
                }
                else {
                        printf("Whoa! Pleading ignorace of CID file type '%c'\n", gParam.CIDfileType);
                        goto problem;
                }
            }
        
/* Read the data */
            
            massToAdd.mOverZ    = -1;    /*test to see if real data entered later*/
            massToAdd.intensity = -1;
            if (gParam.CIDfileType == 'F') /* Finnigan's ICIS text format */
            {
                sscanf(stringBuffer, "%*d %f %d", &massToAdd.mOverZ, &massToAdd.intensity);
            }    
            else if (gParam.CIDfileType == 'T') /* Tab text format */
            {
                sscanf(stringBuffer, "%f %d", &massToAdd.mOverZ, &massToAdd.intensity);
                if (i == 1 && (massToAdd.mOverZ <= 1 || massToAdd.mOverZ > 1000))
                {
                    printf("The datafile does not appear to be in tab-delimited (T) format\n"
                           "Be sure the CID file type is set correctly in the Lutefisk.params file.\n");
                    goto problem;
                }
            }
            else if (gParam.CIDfileType == 'Q') /* Micromass' pkl '.dta' text format */
            {
                sscanf(stringBuffer, "%f %f", &massToAdd.mOverZ, &intensityAsReal);
                massToAdd.intensity = (INT_4) intensityAsReal;
                if (i == 1 && (massToAdd.mOverZ <= 1 || massToAdd.mOverZ > 1000))
                {
                    printf("The datafile does not appear to be in Micromass' pkl '.dta' (Q) format\n"
                           "Be sure the CID file type is set correctly in the Lutefisk.params file.\n");
                    goto problem;
                }
            }
            else if (gParam.CIDfileType == 'L') /* Finnigan's LCQ text format */
            {
                if (stringBuffer[0] == '\n')
                {
                    continue;    /*skip any blank lines*/
                }
                sscanf(stringBuffer, "%s", stringBuffer2);
                if (!strcmp(stringBuffer2, "saturated"))
                {
                    continue;    /*finnigan sticks this into every other line for some reason*/
                }
                for (j = 0; j < 258; j++)
                {
                    if(stringBuffer[j] == ',')    /*replaces commas with spaces*/
                    {
                        stringBuffer[j] = ' ';
                    }
                }
                sscanf(stringBuffer, "%*[Packet] %*[#] %*d %*[intensity] %*[=] %f %*[mass/position] %*[=] %f", &intensityAsReal, &massToAdd.mOverZ);
                massToAdd.intensity = (INT_4) intensityAsReal;
            }
            else if (gParam.CIDfileType == 'D') /* Finnigan's '.dta' text format */
            {
                                /* In case the intensity is a real value, we will read it this way. */
                sscanf(stringBuffer, "%f %f", &massToAdd.mOverZ, &intensityAsReal);
                massToAdd.intensity = (INT_4) intensityAsReal;
            }
  
            
            if (massToAdd.mOverZ < -1 || massToAdd.intensity < 0)
            {
                printf("There is something wrong with the data file.\n");
                goto problem;
            }
            
            if (firstNumberFlag)
            {
                msms.scanMassLow = massToAdd.mOverZ;
                firstNumberFlag = false;
            }
            if(massToAdd.mOverZ > msms.scanMassHigh)
            	msms.scanMassHigh = massToAdd.mOverZ;
                
            if (oldMassValue == massToAdd.mOverZ)
            {
                if (massToAdd.intensity >= oldIonIntensity)
                {
                    MSDataList->mass[MSDataList->numObjects - 1] = massToAdd;
                    oldMassValue = massToAdd.mOverZ;
                    oldIonIntensity = massToAdd.intensity;
                }
            }
            else
            {
                if (!AddToList(&massToAdd, MSDataList)) 
                {
                    printf("Ran out of room for datapoints!\n");
                    goto problem;
                }
    
                oldMassValue = massToAdd.mOverZ;
                oldIonIntensity = massToAdd.intensity;
            }
        }

        fclose(fp);
    }
 
    free(stringBuffer);
    free(stringBuffer2);
 
    return MSDataList;
  
  problem:
  
        printf("Quitting.");
        exit(1);  
        return(NULL);
}
    
  
/************************************ CalibrationCorrection ********************************
*
*    If sufficiently intense immonium ions and y1 ions (for tryptic peptides) are available,
*    determine a mass correction to be applied to the list of ions.
*
*/
void CalibrationCorrection(tMSDataList *inPeakList)
{
    REAL_4 immoniumIons[13], y1Arg, y1Lys, intensityCutoff, calibrationMass[15], calibrationIntensity[15];
    REAL_4 offsetMass, totalOffsetIntensity, lysErr, argErr;
    INT_4  i, ionNum;
    
    tMSData *currPtr = NULL;
    tMSData *ptrOfNoReturn = NULL;
    
    /*Calculate the values for y1Arg, y1Lys, and immoniumIons.*/
    y1Arg = gMonoMass[R] + 3 * gElementMass[HYDROGEN] + gElementMass[OXYGEN];
    y1Lys = gMonoMass[K] + 3 * gElementMass[HYDROGEN] + gElementMass[OXYGEN];
    immoniumIons[0] = gMonoMass[P] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[1] = gMonoMass[V] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[2] = gMonoMass[L] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[3] = gMonoMass[M] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[4] = gMonoMass[H] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[5] = gMonoMass[F] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[6] = gMonoMass[Y] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[7] = gMonoMass[W] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[8] = gMonoMass[T] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[9] = gMonoMass[S] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[10] = gMonoMass[N] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[11] = gMonoMass[D] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    immoniumIons[12] = gMonoMass[E] + gElementMass[HYDROGEN] - gElementMass[OXYGEN] - gElementMass[CARBON];
    
    
    if (inPeakList->numObjects == 0) return;
    
    /*Calculate the intensity cutoff.*/
    intensityCutoff = 0;
    currPtr = &inPeakList->mass[0];
    ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        intensityCutoff += currPtr->intensity;
        currPtr++;
    }
    
    intensityCutoff = intensityCutoff / inPeakList->numObjects;    /*calc the average intensity*/
    intensityCutoff = intensityCutoff / 4;    /*Only use ions that are greater than 1/4 of the average intensity*/
    
    /*Check if any of the immonium ions or tryptic y1 ions are present.*/
    
    ionNum = 0;    /*reset the ion counter*/
    if(gParam.proteolysis == 'T')    /*look for y1 ions for Arg and Lys if a tryptic peptide.*/
    {
        currPtr = &inPeakList->mass[0];
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->mOverZ > y1Arg + gParam.fragmentErr)
            {
                break;
            }
            if(currPtr->mOverZ <= y1Lys + gParam.fragmentErr && currPtr->mOverZ >= y1Lys - gParam.fragmentErr
                && currPtr->intensity > intensityCutoff)
            {
                calibrationMass[ionNum] = y1Lys - currPtr->mOverZ;
                calibrationIntensity[ionNum] = currPtr->intensity;
                ionNum++;
            }
            if(currPtr->mOverZ <= y1Arg + gParam.fragmentErr && currPtr->mOverZ >= y1Arg - gParam.fragmentErr
                && currPtr->intensity > intensityCutoff)
            {
                calibrationMass[ionNum] = y1Arg - currPtr->mOverZ;
                calibrationIntensity[ionNum] = currPtr->intensity;
                ionNum++;
            }
            currPtr++;
        }
        if(ionNum == 2)    /*if both K and R y1 ions found, then choose the one with the least error*/
        {
            lysErr = calibrationMass[0] - y1Lys;
            argErr = calibrationMass[1] - y1Arg;
            if(lysErr < 0)
            {
                lysErr = lysErr * -1;
            }
            if(argErr < 0)
            {
                argErr = argErr * -1;
            }
            if(lysErr < argErr)
            {
                ionNum = 1;
            }
            else
            {
                calibrationMass[0]      = calibrationMass[1];
                calibrationIntensity[0] = calibrationIntensity[1];
                ionNum = 1;
            }
        }
    }
    if(gParam.proteolysis == 'K')    /*Look for y1 ion of Lys if a Lys-C peptide.*/
    {
        currPtr = &inPeakList->mass[0];
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->mOverZ > y1Lys + gParam.fragmentErr)
            {
                break;
            }
            if(currPtr->mOverZ <= y1Lys + gParam.fragmentErr && currPtr->mOverZ >= y1Lys - gParam.fragmentErr
                && currPtr->intensity > intensityCutoff)
            {
                calibrationMass[ionNum]      = y1Lys - currPtr->mOverZ;
                calibrationIntensity[ionNum] = currPtr->intensity;
                ionNum++;
            }
            currPtr++;
        }
    }
    
    /*Now look for the immonium ions.*/
    for(i = 0; i < 13; i++)
    {
        currPtr = &inPeakList->mass[0];
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->mOverZ > immoniumIons[i] + gParam.fragmentErr)
            {
                break;
            }
            if(currPtr->mOverZ >= immoniumIons[i] - gParam.fragmentErr &&
                currPtr->mOverZ <= immoniumIons[i] + gParam.fragmentErr
                && currPtr->intensity > intensityCutoff)
            {
                calibrationMass[ionNum]      = immoniumIons[i] - currPtr->mOverZ;
                calibrationIntensity[ionNum] = currPtr->intensity;
                ionNum++;
            }
            currPtr++;
        }
    }
    
    /*Calculate the offset values*/ 
    
    offsetMass = 0;
    totalOffsetIntensity = 0;
    for(i = 0; i < ionNum; i++)
    {
        offsetMass += calibrationMass[i] * calibrationIntensity[i];
        totalOffsetIntensity += calibrationIntensity[i];
    }
    if(totalOffsetIntensity == 0)
        return;    /*nothing was found to adjust calibration, so return w/o modifying the masses.*/
    offsetMass = offsetMass / totalOffsetIntensity;    /*obtain the average offset weighted for intensity*/
    offsetMass = offsetMass * 0.5;    /*adjust by only half as much as calculated (dont be too radical)*/
    
    /*Apply the calculated offset value to all of the ion masses.*/
    
    currPtr = &inPeakList->mass[0];
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->mOverZ = currPtr->mOverZ + offsetMass;
        currPtr++;
    }
    
    msms.scanMassLow = msms.scanMassLow + offsetMass;
    msms.scanMassHigh = msms.scanMassHigh + offsetMass;
    
    printf("The QTof calibration offset is %f \n", offsetMass);

    return;
}

/***********************************CheckConnections****************************************
*
*    This is where I make sure that ions can be connected to other ions via single amino 
*    acid jumps.  The mass accuracy required is increased for this determination.
*/
void CheckConnections(tMSDataList *inPeakList)
{
    tMSData *currPtr;
    tMSData *ptrOfNoReturn;
    tMSData *nextPtr;
    REAL_4 massDiff, error;
    REAL_4 highMassYIon, lowMassYIon, highMassBIon;
    REAL_4 precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) 
                        / gParam.chargeState;
    INT_4 j, i, maxCharge;
    char test;
    
    if(gParam.fragmentErr >= 0.75 || (gParam.chargeState > 2 && gParam.maxent3 == FALSE))
    {
        return;    /*Invoke this routine only by setting the fragment error to less than 0.75.*/
    }
    
    error = gParam.fragmentErr * 0.5;   /*9/23/03 changed rsj*/

    /* make all indexes negative, later any ions that can connect are made positive, 
       and in the end the ions that continue to have negative indexes are removed */
    currPtr = &inPeakList->mass[0];    
    ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];    
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->index = -1;
        currPtr++;
    }
    
/*    For doubly charged precursors from ion trap data, don't eliminate fragment ions that could be 
    doubly-charged.*/
    if(gParam.fragmentPattern == 'L' && gParam.chargeState == 2)
    {
        currPtr = &inPeakList->mass[0];    
        ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];    
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->mOverZ >= precursor - gMonoMass[W] - gParam.fragmentErr && 
                currPtr->mOverZ <= precursor - 2 * gParam.fragmentErr)
            {
                currPtr->index = 1;
            }
            currPtr++;
        }
    }
    
    /*now find connecting ions*/
    currPtr = &inPeakList->mass[0];    
    while(currPtr < ptrOfNoReturn - 1)
    {
        nextPtr = currPtr + 1;
        while(nextPtr < ptrOfNoReturn)
        {
            massDiff = nextPtr->mOverZ - currPtr->mOverZ;
            if(massDiff > gMonoMass[W] + error)
            {
                break;    /*stop looking if the difference is more than Trp*/
            }
            test = FALSE;
            for(i = 0; i < gAminoAcidNumber; i++)
            {
                if(massDiff <= gMonoMass[i] + error && 
                    massDiff >= gMonoMass[i] - error)
                {
                    test = TRUE;
                    break;
                }
            }
            if(test)
            {
                if(nextPtr->index == -1)
                {
                    nextPtr->index = 1;
                }
                if(currPtr->index == -1)
                {
                    currPtr->index = 1;
                }
            }
            nextPtr++;
        }
        currPtr++;
    }
        
    /*Keep ions that could be terminal y ions or high mass b ions.*/
    
    if(gParam.maxent3)
    {
        maxCharge = 1;  /*maxCharge is one for maxent3 data, since all ions converted to +1*/
    }
    else
    {
        maxCharge = gParam.chargeState;
    }
    
    currPtr = &inPeakList->mass[0];    
    while(currPtr < ptrOfNoReturn)
    {
        for(i = 1; i <= maxCharge; i++)
        {
            for(j = 0; j < gAminoAcidNumber; j++)
            {
                /*Calculate the high mass y ion.*/
                highMassYIon = gParam.peptideMW - gMonoMass[j] + gParam.modifiedNTerm - gElementMass[HYDROGEN];
                highMassYIon = (highMassYIon + (gElementMass[HYDROGEN] * i)) / i;
                
                /*Calculate the low mass y ion.*/
                lowMassYIon = gMonoMass[j] + gElementMass[HYDROGEN] + gParam.modifiedCTerm;
                lowMassYIon = (lowMassYIon + (gElementMass[HYDROGEN] * i)) / i;
                
                /*Calculate the high mass b ion.*/
                highMassBIon =  gParam.peptideMW - gMonoMass[j] - gParam.modifiedCTerm;
                highMassBIon = (highMassBIon + (gElementMass[HYDROGEN] * (i - 1))) / i;
                
                if(currPtr->mOverZ >= highMassYIon - gParam.fragmentErr &&
                    currPtr->mOverZ <= highMassYIon + gParam.fragmentErr)
                {
                    if(currPtr->index == -1)
                    {
                        currPtr->index = 1;    /**Fixed by RSJ, JAT had currptr->intensity*/
                    }
                }
                if(currPtr->mOverZ >= lowMassYIon - gParam.fragmentErr &&
                    currPtr->mOverZ <= lowMassYIon + gParam.fragmentErr)
                {
                    if(currPtr->index == -1)
                    {
                        currPtr->index = 1;    /**Fixed by RSJ, JAT had currptr->intensity*/
                    }
                }
                if(currPtr->mOverZ >= highMassBIon - gParam.fragmentErr &&
                    currPtr->mOverZ <= highMassBIon + gParam.fragmentErr)
                {
                    if(currPtr->index == -1)
                    {
                        currPtr->index = 1;    /**Fixed by RSJ, JAT had currptr->intensity*/
                    }
                }
            }
        }
    
        currPtr++;
    }
    
    
    /*Keep the low mass ions.*/
    currPtr = &inPeakList->mass[0];    
    while(currPtr < ptrOfNoReturn)
    {
        if(currPtr->mOverZ < 148 || (currPtr->mOverZ > 158.5 && currPtr->mOverZ < 159.5))
        {
            if(currPtr->index == -1)
            {
                currPtr->index = 1;
            }
        }
        currPtr++;
    }
    
    /*get rid of the peaks with neg indexes*/
    currPtr = &inPeakList->mass[0];    
    while(currPtr < ptrOfNoReturn)
    {
        if(currPtr->index == -1)
        {
            RemoveFromList(currPtr - &inPeakList->mass[0], inPeakList);
            ptrOfNoReturn--;
            currPtr--;
        }

        currPtr++;
    }            


    /* Reset index values */    
    for (i = 0;    i <    inPeakList->numObjects; i++)
    {
        inPeakList->mass[i].index = i;
    }    
    
    return;
}

/***********************************CheckSignalToNoise**************************************
*
*    This function verifies that the remaining ions selected as "interesting" have an adequate
*    signal to noise ratio.  This ratio is #defined in LutefiskDefinitions.  The noise is
*    calculated as the average intensity 50 u above and 50 u below the ion.  For simplicity,
*    the first element in the linked list is never deleted.
*/

void CheckSignalToNoise(tMSDataList *inPeakList, tMSDataList *inMSDataList)
{
    tMSDataList *neighborhoodList = NULL;
    tMSData *currPtr = NULL;
    tMSData    *previousPtr = NULL;
    tMSData *ptrOfNoReturn = NULL;
    tMSData *currPeakPtr = NULL;
    tMSData *peakPtrOfNoReturn = NULL;
    INT_4     dataNum;
    INT_4     range = 50;    /*This is the range (+/- this num) over which the noise is determined*/
    INT_4     deltaMassNum;
    INT_4    calcDataNum;
    REAL_4     signalToNoise, intensity, noise, deltaMass, deltaMassSum, firstDataPoint;

    
    
/*    
    I ran in to trouble by assuming that profile data would have a continuous non-zero intensity.
    However, it seems that Sciex data and probably others has long stretches of zeroed data.
    This had the effect of having a higher calculated noise value, which caused real ions
    to be excluded using the original s/n exclusion criteria established for continuous
    non-zero data.  To overcome this, I need to figure out what the average data spacing is, 
    and then use this to figure out the noise level.
*/
    
    if (inMSDataList->numObjects < 2) return;
    if (inPeakList->numObjects < 2) return;

    previousPtr = &inMSDataList->mass[0];
    currPtr = &inMSDataList->mass[1];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    deltaMass = 0;
    deltaMassSum = 0;
    deltaMassNum = 0;
    
    
    if(gParam.CIDfileType == 'X' || gParam.CIDfileType == 'D')    /*QTof .dta data*/
    {
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->intensity != 0 && previousPtr->intensity != 0)    /*only look at non-zero data points*/
            {
                deltaMass = currPtr->mOverZ - previousPtr->mOverZ;
                if(deltaMass <= 2)    /*dta files seem to have "ions" spaced every u or so*/
                {
                    deltaMassSum += deltaMass;
                    deltaMassNum++;
                }
            }
            currPtr++;
            previousPtr++;
        }
        
        if(deltaMass <= 0 || deltaMassNum < 25)    
        {
            return;    /*if the dta data has been created so that there are few ions less than 2 da apart,
                    then its difficult to assess the noise level*/
        }
        deltaMass = deltaMassSum / (REAL_4)deltaMassNum;
    }
    else
    {
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->intensity != 0 && previousPtr->intensity != 0)    /*only look at non-zero data points*/
            {
                deltaMass = currPtr->mOverZ - previousPtr->mOverZ;
                if(deltaMass <= (gParam.peakWidth * 0.5))
                {
                    deltaMassSum += deltaMass;
                    deltaMassNum++;
                }
            }
            currPtr++;
            previousPtr++;
        }

        if(deltaMassNum <= 0)    /*if deltaMass not found try setting the window a bit wider*/
        {
            previousPtr = &inMSDataList->mass[0];
            currPtr = &inMSDataList->mass[1];
            deltaMass = 0;
            deltaMassSum = 0;
            deltaMassNum = 0;
            while(currPtr < ptrOfNoReturn)
            {
                if(currPtr->intensity != 0 && previousPtr->intensity != 0)    /*only look at non-zero data points*/
                {
                    deltaMass = currPtr->mOverZ - previousPtr->mOverZ;
                    if(deltaMass <= gParam.peakWidth)
                    {
                        deltaMassSum += deltaMass;
                        deltaMassNum++;
                    }
                }
                currPtr++;
                previousPtr++;
            }
        }
        if(deltaMassNum == 0)
            return;    /*avoid divide by zero*/
        deltaMass = deltaMassSum / (REAL_4)deltaMassNum;
    }

    if(deltaMass <= 0)
        return;    /*to prevent divide by zero below*/
    
    
    neighborhoodList = (tMSDataList *) CreateNewList( sizeof(tMSData), 1000, 1000 );
    if (!neighborhoodList) 
    {
        printf("Ran out of memory in CheckSignalToNoise()!\n");
        exit(1);
    }


/*    Start looking at each peak*/
    currPeakPtr = &inPeakList->mass[0];
    peakPtrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];
    while(currPeakPtr < peakPtrOfNoReturn)
    {
/*    First find the ion's neighbors, excluding the ion itself*/
        firstDataPoint = 0;
        dataNum = 0;
        calcDataNum = 0;
        intensity = 0;
        signalToNoise = 0;
        currPtr = &inMSDataList->mass[0];
        
        while(currPtr < ptrOfNoReturn)
        {
            if(currPtr->intensity != 0)    /*dont look at data points w/ zero intensity*/
            {
                if(currPtr->mOverZ > currPeakPtr->mOverZ + range)
                {
                    break;
                }
                else 
                {
                    if (currPtr->mOverZ > currPeakPtr->mOverZ - range
                        && (currPtr->mOverZ < currPeakPtr->mOverZ - gParam.peakWidth
                            || currPtr->mOverZ > currPeakPtr->mOverZ + gParam.peakWidth))
                    {
                        if (!AddToList(currPtr, neighborhoodList))
                        {
                            printf("Ran out of memory in CheckSignalToNoise()!\n");
                            exit(1);
                        }
                    }
                }
            }
            currPtr++;
        }
        
    
        if(neighborhoodList->numObjects > 0)
        {    
            /* Sort the neighborhoodList in order of decreasing intensity */
            qsort(neighborhoodList->mass,(size_t)neighborhoodList->numObjects,
                  (size_t)sizeof(tMSData),IntensityDescendSortFunc);
        
            /* The noise is the median intensity in the list */    
            noise = neighborhoodList->mass[(INT_4)(neighborhoodList->numObjects/2)].intensity;
            
            signalToNoise = (currPeakPtr->intensity) / noise;
        }
        else
        {
            signalToNoise = 100;    /*if there are no data points for measuring noise then give high s/n*/
        }
        
        /* less than 50 data points used to determine the noise level is considered insufficient */
        if ((signalToNoise < SIGNAL_NOISE && neighborhoodList->numObjects > 50) 
            || signalToNoise == 0)    
        {
            /*dont get rid of immonium ions*/
            if(currPeakPtr->mOverZ > 148 && (currPeakPtr->mOverZ < 158.5 || currPeakPtr->mOverZ > 159.5))
            {
                RemoveFromList(currPeakPtr - &inPeakList->mass[0], inPeakList);
                peakPtrOfNoReturn--;
                currPeakPtr--;
            }
        }
    
        neighborhoodList->numObjects = 0;
        currPeakPtr++;
    }
    
    if (neighborhoodList) DisposeList(neighborhoodList);
    
    return;
}

/***********************************GuessAtTheFragmentPattern*******************************
*
*    If the default ("D") is used for the fragment pattern, then the program takes a stab at
*    figuring out if the data is from an LCQ or a TSQ.  It does this by assuming that TSQ
*    data always starts at a mass lower than 0.2 x the precursor m/z.  If changes occur to 
*    the LCQ that permit lower start masses in MS/MS data, then this will need to be changed.
*
*/
void GuessAtTheFragmentPattern()
{
    REAL_4 precursor;
    
    precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / 
                gParam.chargeState;
                
    precursor = precursor * 0.2;
    
    if(precursor > msms.scanMassLow)
    {
        if(gParam.fragmentErr >= 0.15 * gMultiplier)
        {
            gParam.fragmentPattern = 'T'; /* TSQ (triple quad) */
        }
        else
        {
            gParam.fragmentPattern = 'Q'; /* Q-TOF */
        }
    }
    else
    {
        gParam.fragmentPattern = 'L'; /* LCQ (ion trap) */
    }


    return;
}

/***********************************DefectCorrection****************************************
*
*    The observed mass defect is compared to the theoretical mass defect, and corrections are
*    made.  I've found that my LCQ data can have slightly lower than expected masses at higher
*    m/z; for example, an ion at 1100.2 is really at 1100.6.  Since the low m/z end is often
*    ok, I cannot use a mass offset across the entire m/z range.  This more intelligent defect
*    correction should allow for the use of tighter error tolerances (+/- 0.5 Da).
*
*/
void DefectCorrection(tMSDataList *inPeakList)
{
    tMSData     *currPtr = NULL;
    tMSData     *previousPtr = NULL;
    tMSData     *ptrOfNoReturn = NULL;
    INT_4 integerMass[200], ionNum, i;
    REAL_4 mass[200], defect[200], precursor; 
    REAL_8 aObserved, bObserved, sumOfXSquared, sumOfY, sumOfX;
    REAL_8 sumOfXTimesY, aTheory = 0, bTheory = 0.00050275;
    REAL_8 numerator, denominator, observedDefect, theoryDefect, additionalDefect, testMass;
    char doublyCharged[200];
    
    precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / 
                gParam.chargeState;
    
/*    Create mass and integerMass arrays from linked list data*/
    ionNum = 0;
    currPtr = &inPeakList->mass[0];
    ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        mass[ionNum] = currPtr->mOverZ;
        theoryDefect = bTheory * mass[ionNum];          /*added 8/5/03*/
        integerMass[ionNum] = mass[ionNum] - theoryDefect + 0.5;        /*added -theoryDefect + 0.5 8/5/03*/
        defect[ionNum] = mass[ionNum] - integerMass[ionNum];
        ionNum++;
        if(ionNum >= 200)
            return;    /*if too many ions, then return w/o making corrections*/
        currPtr++;
    }
    
/*    Find potential doubly-charged ions, and mark them so that they are not changed.*/
    for(i = 0; i < ionNum; i++)
    {
        doublyCharged[i] = 0;    /*initialize*/
    }
    if(gParam.chargeState == 2)    
    {
        for(i = 0; i < ionNum; i++)
        {
            if(mass[i] > precursor - (3 * gMonoMass[W] / 2) - gParam.fragmentErr &&
                mass[i] < precursor + gParam.fragmentErr)
            {
                doublyCharged[i] = 1;
            }
        }
            /*for(j = 0; j < ionNum; j++)
            {
                if(((mass[i] * 2) - gElementMass[HYDROGEN] <= mass[j] + gParam.fragmentErr)
                    && ((mass[i] * 2) - gElementMass[HYDROGEN] >= mass[j] - gParam.fragmentErr))
                {
                    doublyCharged[i] = 1;
                }
            }*/
    }
/*    Calculated the mass defect for each data point*/
   /* for(i = 0; i < ionNum; i++)
    {
        defect[i] = mass[i] - integerMass[i];
        if(integerMass[i] < 700 && doublyCharged[i] == 0)    
        {
            if(defect[i] > 0.7)    //for low mass, if error is in other direction, the defect should be negative.
            {
                defect[i] = 0;
                integerMass[i] = mass[i] + 0.5;
                mass[i] = integerMass[i];    //bump the value up to at least the integer value
            }
        }
    }*/
    
/*    Now for the least squares calculation of a straight line.*/
    sumOfX = 0;    /*calc sumOfX*/
    for(i = 0; i < ionNum; i++)
    {
        sumOfX += mass[i];
    }
    
    sumOfY = 0;    /*calc sumOfY*/
    for(i = 0; i < ionNum; i++)
    {
        sumOfY += defect[i];
    }
    
    sumOfXSquared = 0;    /*calc sumOfXSquared*/
    for(i = 0; i < ionNum; i++)
    {
        sumOfXSquared = sumOfXSquared + (mass[i] * mass[i]);
    }
    
    sumOfXTimesY = 0;    /*calc of sumOfXTimesY*/
    for(i = 0; i < ionNum; i++)
    {
        sumOfXTimesY = sumOfXTimesY + (mass[i] * defect[i]);
    }
    
/*    From equation y = bx + a, calc a first*/
    numerator = (sumOfXSquared * sumOfY) - (sumOfX * sumOfXTimesY);
    denominator = (ionNum * sumOfXSquared) - (sumOfX * sumOfX);
    if(denominator == 0)
    {
        printf("DefectCorrection:  denominator = 0! Quitting.\n");
        exit(1);
    }
    aObserved = numerator / denominator;
    
/*    Now calculate b*/
    numerator = (ionNum * sumOfXTimesY) - (sumOfX * sumOfY);
    denominator = (ionNum * sumOfXSquared) - (sumOfX * sumOfX);
    if(denominator == 0)
    {
        printf("DefectCorrection:  denominator = 0\n");
        exit(1);
    }
    bObserved = numerator / denominator;
    
/*    Now make the mass corrections.*/
    for(i = 0; i < ionNum; i++)
    {
        if(mass[i] > 0 && doublyCharged[i] == 0)
        {
            testMass = (mass[i] * bTheory) - defect[i];
            if(testMass > 0.15)    /*if already close, then don't bother*/
            { 
                observedDefect = (bObserved * mass[i]) + aObserved;
                theoryDefect = (bTheory * mass[i]) + aTheory;
                additionalDefect = theoryDefect - observedDefect;
                if(additionalDefect < 0)
                {
                    additionalDefect = (theoryDefect - defect[i]) / 2;
                }
                if(additionalDefect > 0)
                {
                    mass[i] = mass[i] + additionalDefect;       
                }
            }
            testMass = defect[i] - (mass[i] * bTheory);
            if(testMass > 0.15) 
            {
                observedDefect = (bObserved * mass[i]) + aObserved;
                theoryDefect = (bTheory * mass[i]) + aTheory;
                additionalDefect = theoryDefect - observedDefect;
                if(additionalDefect > 0)
                {
                    additionalDefect = (theoryDefect - defect[i]) / 2;
                }
                if(additionalDefect < 0)
                {
                    mass[i] = mass[i] + additionalDefect;
                }
            }
        }
    }
    
    currPtr = &inPeakList->mass[0];
    i = 0;
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->mOverZ = mass[i];
        i++;
        currPtr++;
    }


    return;
}

/***********************************NormalizeIntensity**************************************
*
*    This function normalizes the CID data intensity to the fourth most intense ion.  It
*    seems to not be unusual for LCQ data to have a few favored fragmentation pathways that
*    result in a couple of intense ions.  In order to increase the spread in the scoring of
*    candidate sequences, these intense ions are reduced to the intensity of the fourth most
*    abundant ion.
*
*/
void NormalizeIntensity(tMSDataList *inMSDataList)
{

    /* Sort the MSData in order of decreasing intensity */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),IntensityDescendSortFunc);

    inMSDataList->mass[0].intensity = inMSDataList->mass[3].intensity;
    inMSDataList->mass[1].intensity = inMSDataList->mass[3].intensity;
    inMSDataList->mass[2].intensity = inMSDataList->mass[3].intensity;
    
    /* Resort the MSData in order of increasing mass */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),MassAscendSortFunc);
    

    return;
}

/***********************************CentroidOrProfile***************************************
*
*    The mass differences between adjacent data points are calculated, and a standard deviation
*    for these differences is used to establish if the data is profile (very little deviation)
*    or centroided data (large deviation).
*
*/
void CentroidOrProfile(tMSDataList *inMSDataList)
{
    tMSData     *currPtr;
    tMSData     *prevPtr;
    tMSData     *ptrOfNoReturn;
    INT_4         threshold = 0;
    INT_4        numberOfMassDiffs = 0;
    INT_4        i;
    REAL_8         massDiff[100];
    REAL_8        massDiffAv = 0;
    REAL_8        standardDeviation = 0;
    
    
    if (inMSDataList->numObjects < 2) return; 
        
    threshold = FindThreshold(inMSDataList);

    prevPtr = &inMSDataList->mass[0];
    currPtr = &inMSDataList->mass[1];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
        
    while (currPtr < ptrOfNoReturn && numberOfMassDiffs < 100) 
    {

        if(currPtr->intensity > threshold && prevPtr->intensity > threshold)
        {
            massDiff[numberOfMassDiffs] = (currPtr->mOverZ) - (prevPtr->mOverZ);
            massDiffAv += massDiff[numberOfMassDiffs];
            numberOfMassDiffs++;  
        }
        currPtr++;
        prevPtr++;
    }
    
    if(numberOfMassDiffs == 0) return; /* Avoid potential divide-by-zero */
    
    massDiffAv = massDiffAv / numberOfMassDiffs;
    
    for(i = 0; i < numberOfMassDiffs; i++)
    {
        standardDeviation = standardDeviation + 
                            ((massDiff[i] - massDiffAv) * (massDiff[i] - massDiffAv));
    }
    standardDeviation = standardDeviation / numberOfMassDiffs;
    standardDeviation = sqrt(standardDeviation);
    
    if(massDiffAv < 1)
    {
        if(standardDeviation < 0.5)
        {
            gParam.centroidOrProfile = 'P';
        }
        else
        {
            gParam.centroidOrProfile = 'C';
        }
    }
    else
    {
        gParam.centroidOrProfile = 'C';
    }

    return;
}

/***********************************FindBYGoldenBoys****************************************
*
*       For LCQ data, find b/y ion pair complements that represent cleavage at the same amide
*       bond.  Designate these as golden boy ions that are difficult to get rid of simply on the
*       basis of ion intensity.
*
*/

void FindBYGoldenBoys(tMSDataList *inMSDataList)
{
        REAL_4  *massList, *mass2List, testMass, *pairMass, avePairMass;
        REAL_8  stDev;
        INT_4   maxIonNum = gGraphLength / gMultiplier;
        INT_4   ionNum, i, j, pairNum;
        char    *goodOrBad;
        
        tMSData         *currPtr = NULL;
    tMSData             *ptrOfNoReturn = NULL;
    
/*    Set aside some space for these arrays.*/
    massList = (float *) malloc(maxIonNum * sizeof(REAL_4));
    if(massList == NULL)
    {
        printf("FindBYGoldenBoys:  Out of memory.");
        exit(1);
    }
    mass2List = (float *) malloc(maxIonNum * sizeof(REAL_4));
    if(mass2List == NULL)
    {
        printf("FindBYGoldenBoys:  Out of memory.");
        exit(1);
    }
    pairMass = (float *) malloc(maxIonNum * sizeof(REAL_4));
    if(massList == NULL)
    {
        printf("FindBYGoldenBoys:  Out of memory.");
        exit(1);
    }
    goodOrBad = (char *) malloc(maxIonNum * sizeof(char));
    if(goodOrBad == NULL)
    {
        printf("FindBYGoldenBoys:  Out of memory.");
        exit(1);
    }
    
/*    Initialize some variables.*/
        currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->normIntensity = 0;        /*set normIntensity field to zero*/
        currPtr++;
    }
    
    ionNum              = 0;
    pairNum     = 0;
    avePairMass = 0;
    stDev               = 0;
    
    for(i = 0; i < maxIonNum; i++)
    {
        massList[i]  = 0;
        mass2List[i] = 0;
        goodOrBad[i] = 0;
    }
        
/*    Fill in the mass array assuming singly-charged ions.*/
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        massList[ionNum] = currPtr->mOverZ;
        ionNum++;
        currPtr++;
    }

/*    Fill in the mass array assuming doubly-charged ions.*/
        if(gParam.chargeState > 2)
        {
                for(i = 0; i < ionNum; i++)
                {
                        testMass = massList[i] * 2 - gElementMass[HYDROGEN];
                        if(testMass < gParam.peptideMW - gMonoMass[G] + gParam.fragmentErr &&
                                testMass > 700) /*doubly charged ions have to be in the right mass range*/
                        {
                        mass2List[i] = testMass;
                    }
                    else
                    {
                        mass2List[i] = 0;       /*zero is a flag that it could not be a doubly-charged ion*/
                    }
                }
        }

/*      Find a suitable error*/
        /*assume all ions are singly-charged*/
        for(i = 0; i < ionNum - 1; i++)
        {
                for(j = i + 1; j < ionNum; j++)
                {
                        testMass = massList[i] + massList[j] - 2 * gElementMass[HYDROGEN];
                        if(testMass <= gParam.peptideMW + gParam.peptideErr * 2 &&
                                testMass >= gParam.peptideMW - gParam.peptideErr * 2)
                        {
                                pairMass[pairNum] = testMass;   /*collect the data*/
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
                                if(mass2List[j] > massList[i])
                                {
                                        testMass = massList[i] + mass2List[j] - 2 * gElementMass[HYDROGEN];
                                        if(testMass <= gParam.peptideMW + gParam.peptideErr * 2 &&
                                                testMass >= gParam.peptideMW - gParam.peptideErr * 2)
                                        {
                                                pairMass[pairNum] = testMass;   /*collect the data*/
                                                pairNum++;
                                        }
                                }
                        }
                }
        }
        
        /*now calculate the stDev error*/
        if(pairNum < 3)
        {
                stDev = gParam.peptideErr;      /*not enough pairs of ions to determine standard deviation
                                                                        so it gets defined as the peptide error from the params file*/
        }
        else    /*enough data to take a stab at finding standard deviation*/
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
                stDev = 2 * gParam.peptideErr;  /*don't let the error be too big*/
        }
        else if(stDev < 0.5 * gParam.peptideErr)
        {
                stDev = 0.5 * gParam.peptideErr;        /*or too small*/
        }


        /*find pairs of masses that are close to the peptide molecular weight*/
        /*first assume the ions are all singly-charged*/
        pairNum = 0;
        for(i = 0; i < ionNum - 1; i++)
        {
                for(j = i + 1; j < ionNum; j++)
                {
                        testMass = massList[i] + massList[j] - 2 * gElementMass[HYDROGEN];
                        if(testMass <= gParam.peptideMW + stDev &&
                                testMass >= gParam.peptideMW - stDev)
                        {
                                goodOrBad[i] = 1;
                                goodOrBad[j] = 1;
                                pairNum++;
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
                                if(mass2List[j] > massList[i])
                                {
                                        testMass = massList[i] + mass2List[j] - 2 * gElementMass[HYDROGEN];
                                        if(testMass <= gParam.peptideMW + stDev &&
                                                testMass >= gParam.peptideMW - stDev)
                                        {
                                                goodOrBad[i] = 1;
                                                goodOrBad[j] = 1;
                                                pairNum++;
                                        }
                                }
                        }
                }
        }
        
/*
    The normIntensity field for the ms data pointers contain 0 if not a goldenBoy and a 1 if
    it is.
*/

    for(i = 0; i < ionNum; i++)
    {
        if(goodOrBad[i] != 0)
        {
            currPtr = &inMSDataList->mass[0];
            ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
            while(currPtr < ptrOfNoReturn)
            {

                if(massList[i] == currPtr->mOverZ)
                {
                    currPtr->normIntensity = 1;
                    break;
                }
                currPtr++;
            }
        }
    }

/*      free the variables*/
        free(massList);
        free(mass2List);
        free(pairMass);
        free(goodOrBad);
        
        return;
}

/***********************************FindTheGoldenBoys***************************************
*
*    Two arrays are set up, one that contains the ion m/z and another that is either 1 or 0,
*    depending on whether the ion is a golden boy or not (golden boys are ions that can be
*    connected by single amino acid jumps to either 147 or 175, and are m/z less than the 
*    precursor ion).  By searching through the ms data list, I start at 147 and make one aa jumps.
*    If an ion is present then the same index position for the second array is reset to 1.
*
*/
void FindTheGoldenBoys(tMSDataList *inMSDataList)
{
    REAL_4 *massList, precursor, lys, arg, plusArgLys[2 * AMINO_ACID_NUMBER];
    REAL_4 testMass, err;
    char *goodOrBad, test;
    INT_4 i, j, k, ionNum, *intensityList, cutoff, goldenBoyNum;
    INT_4 maxIonNum = gGraphLength / gMultiplier;    /*don't need GRAPH_LENGTH numbers of ions for the
                                                    arrays of massList, goodOrBad, and intensityList*/
    tMSData         *currPtr = NULL;
    tMSData         *ptrOfNoReturn = NULL;
    
    
/*    Set aside some space for these arrays.*/
    massList = (float *) malloc(maxIonNum * sizeof(REAL_4));
    if(massList == NULL)
    {
        printf("FindTheGoldenBoys:  Out of memory.");
        exit(1);
    }
    goodOrBad = (char *) malloc(maxIonNum * sizeof(char));
    if(goodOrBad == NULL)
    {
        printf("FindTheGoldenBoys:  Out of memory.");
        exit(1);
    }
    intensityList = (int *) malloc(maxIonNum * sizeof(INT_4));
    if(intensityList == NULL)
    {
        printf("FindTheGoldenBoys:  Out of memory.");
        exit(1);
    }


/*    Initialize some variables.*/
    lys = gMonoMass[K] + 3 * gElementMass[HYDROGEN] + gElementMass[OXYGEN];
    arg = gMonoMass[R] + 3 * gElementMass[HYDROGEN] + gElementMass[OXYGEN];
    
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->normIntensity = 0;        /*set normIntensity field to zero*/
        currPtr++;
    }
    for(i = 0; i < maxIonNum; i++)
    {
        massList[i] = 0;
        intensityList[i] = 0;
        goodOrBad[i] = 0;
    }
    for(i = 0; i < gAminoAcidNumber; i++)
    {
        plusArgLys[i] = lys + gMonoMass[i];
    }
    j = 0;
    for(i = gAminoAcidNumber; i < 2 * gAminoAcidNumber; i++)
    {
        plusArgLys[i] = arg + gMonoMass[j];
        j++;
    }
    precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN]) / gParam.chargeState;
    ionNum = 0;
    err = gParam.fragmentErr / 4;     /*The error between ions is less than the error of the 
                                    calc vs obsd masses.*/
    
/*    Fill in the mass array.*/
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        if(currPtr->mOverZ > lys - gParam.fragmentErr &&
            currPtr->mOverZ < precursor - gParam.fragmentErr * 2 &&
            currPtr->mOverZ < GOLDEN_BOY_MAX)
        {
            massList[ionNum] = currPtr->mOverZ;
            intensityList[ionNum] = currPtr->intensity;
            ionNum++;
        }
        currPtr++;
    }
    

/*
    Seed the array goodOrBad by putting a 1 in the position that have 147 or 175.  If both
    are absent, then I seed with the array plusArgLys.
*/
    test = 1;
    for(i = 0; i < ionNum; i++)
    {
        if((massList[i] <= lys + gParam.fragmentErr) && 
            (massList[i] >= lys - gParam.fragmentErr))
        {
            goodOrBad[i] = 1;    /*you found 147*/
            test = 0;
            break;
        }
    }
    for(i = 0; i < ionNum; i++)
    {
        if((massList[i] <= arg + gParam.fragmentErr) && 
            (massList[i] >= arg - gParam.fragmentErr))
        {
            goodOrBad[i] = 1;    /*you found 175*/
            test = 0;
            break;
        }
    }
    if(test)
    {
        for(i = 0; i < ionNum; i++)
        {
            for(j = 0; j < 2 * gAminoAcidNumber; j++)
            {
                if((massList[i] <= plusArgLys[j] + gParam.fragmentErr) &&
                    (massList[i] >= plusArgLys[j] - gParam.fragmentErr))
                {
                    goodOrBad[i] = 1;
                }
            }
        }
    }
    
/*
    Start at the low mass end and work up trying to connect y ions.  Anything that can be
    connected to 147 or 175 is given a value of 1 in the goodOrBad array.
*/

    for(i = 0; i < ionNum - 1; i++)
    {
        if(goodOrBad[i] != 0)
        {
            for(j = i + 1; j < ionNum; j++)
            {
                testMass = massList[j] - massList[i];
                if((testMass <= gMonoMass[W] + err) &&
                    (testMass >= gMonoMass[G] - err))
                {
                    for(k = 0; k < gAminoAcidNumber; k++)
                    {
                        if((testMass <= gMonoMass[k] + err) &&
                            (testMass >= gMonoMass[k] - err))
                        {
                                if(massList[i] <= lys + err && massList[i] >= lys - err)
                                {       
                                        goodOrBad[j] = -1;      /*tag y2 ions so that they don't get tossed*/
                                }
                                else if(massList[i] <= arg + err && massList[i] >= arg - err)
                                {
                                        goodOrBad[j] = -1;      /*tag y2 ions so that they don't get tossed*/
                                }
                                else
                                {
                                goodOrBad[j] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

/*      Make sure we don't lose the y1 ions by marking them as -1*/
        for(i = 0; i < ionNum; i++)
        {
                if(massList[i] <= lys + err && massList[i] >= lys - err)
        {       
                goodOrBad[i] = -1;      /*tag y2 ions so that they don't get tossed*/
        }
        if(massList[i] <= arg + err && massList[i] >= arg - err)
        {
                goodOrBad[i] = -1;      /*tag y2 ions so that they don't get tossed*/
        }
    }
                                
/*
    Eliminate low intensity goldenBoys from the goldenBoy list.
*/

    cutoff = 0;
    goldenBoyNum = 0;
    for(i = 0; i < ionNum; i++)
    {
        if(goodOrBad[i] == 1)
        {
            cutoff += intensityList[i];
            goldenBoyNum++;
        }
    }
    if(goldenBoyNum == 0) return;
    cutoff = (cutoff / goldenBoyNum) * GOLDEN_BOY_CUTOFF;
    
    for(i = 0; i < ionNum; i++)
    {
        if(goodOrBad[i] == 1)
        {
            if(intensityList[i] < cutoff)
            {
                 goodOrBad[i] = 0;
            }
        }
    }
    
    /*change the -1 values that mark y2 ions back to +1 values so that they get counted as goldenboys*/
    for(i = 0; i < ionNum; i++)
    {
        if(goodOrBad[i] == -1)
        {
                goodOrBad[i] = 1;
        }
    }
    
/*
    The normIntensity field for the ms data pointers contain 0 if not a goldenBoy and a 1 if
    it is.
*/

    for(i = 0; i < ionNum; i++)
    {
        if(goodOrBad[i] != 0)
        {
            currPtr = &inMSDataList->mass[0];
            ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
            while(currPtr < ptrOfNoReturn)
            {

                if(massList[i] == currPtr->mOverZ)
                {
                    currPtr->normIntensity = 1;
                    break;
                }
                currPtr++;
            }
        }
    }
        
/*    Free the arrays.*/
    free(massList);
    free(goodOrBad);
    free(intensityList);
    return;
}

/********************************GetPeakWidth***********************************************
*
*    GetPeakWidth finds the peak width when the auto-peakWidth option is chosen by setting
*    the peakWidth to zero in the .params file.
*
*/

REAL_4 GetPeakWidth(tMSDataList *inMSDataList)
{

    tMSDataList     *bigTreeList = NULL;
    INT_4            i;
    tMSData         *currPtr;
    tMSData         *ptrOfNoReturn;
    REAL_4            precursor;
    INT_4            topIndex;
    INT_4            halfIntensity;
    REAL_4            slope;
    REAL_4            intercept;
    REAL_4            leadingMass;
    REAL_4            trailingMass;
    REAL_4            peakWidth;
    REAL_4            peakWidthSum = 0;
    REAL_4            peakWidthSquaredSum = 0;
    REAL_4            avgPeakWidth;
    REAL_4            stdDev;
    REAL_4            tolerance;
    
    
    if (inMSDataList->numObjects == 0)
    {
        printf("GetPeakWidth: no data in inMSDataList\n");
        exit(1);
    }
    
    precursor = (gParam.peptideMW + gParam.chargeState) / gParam.chargeState;

    /* Sort the MSData in order of decreasing intensity */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),IntensityDescendSortFunc);


    bigTreeList = (tMSDataList *) CreateNewList( sizeof(tMSData), 10, 1 );
    if (!bigTreeList) 
    {
        printf("Ran out of memory in GetPeakWidth()!\n");
        exit(1);
    }

    /* Make bigTrees the top ten most intense peaks that are not the precursor. */
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    while (currPtr < ptrOfNoReturn && bigTreeList->numObjects < 10) 
    {
        /* Only gather peaks below 600 and don't take the precursor */
        if (currPtr->mOverZ < 600 
            && (currPtr->mOverZ < precursor - 2.0
                || currPtr->mOverZ > precursor + 2.0))
/* XXXXXX WHICH SHOULD I USE????? JAT
//            && (currPtr->mOverZ < precursor - gParam.peptideErr
//                || currPtr->mOverZ > precursor + gParam.peptideErr))
*/
        {
            /* Don't take a peak that overlaps one already on the list */
            for (i = 0; i < bigTreeList->numObjects; i++)
            {
                if (currPtr->mOverZ > bigTreeList->mass[i].mOverZ - 5.0
                    &&     currPtr->mOverZ < bigTreeList->mass[i].mOverZ + 5.0)    
                {
                    break; /* Too close to a tree we already have */
                }
            }
            if (i == bigTreeList->numObjects) 
            {    
                if(!AddToList(currPtr, bigTreeList)) 
                {
                    printf("Ran out of memory in GetPeakWidth()!\n");
                    exit(1);
                }
            }
        }            
        currPtr++;
    }
    
    /* Resort the MSData in order of increasing mass */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),MassAscendSortFunc);


    /* Remove big tree peaks if they are < 10% of the highest peak. */
    for (i = 1; i < bigTreeList->numObjects; i++)
    {
        if (bigTreeList->mass[i].intensity < 0.10 * bigTreeList->mass[0].intensity)
        {
            RemoveFromList(i, bigTreeList);
            i--;
        }
    }
    
    /* Now find the half-height peak width of each remaining bigTree. */
    for (i = 0; i < bigTreeList->numObjects; i++)
    {    
        currPtr = &inMSDataList->mass[0];
        ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
        
        while (currPtr < ptrOfNoReturn) 
        {
            if (currPtr->mOverZ == bigTreeList->mass[i].mOverZ) break;        
            currPtr++;
        }    
        
        topIndex = currPtr - &inMSDataList->mass[0];
        halfIntensity = (INT_4)(0.5 * bigTreeList->mass[i].intensity);

        /*
        //                   peak top
        //                      |
        //     leading          .        trailing
        //       edge          . .         edge
        //                    .   .
        //    ______________ . ___ . ______________ 50% peak height
        //                  .       .
        //                 .         .
        //               .             .
        //     . . . . .     |-----|     . . . . . 
        //                  peak width
        */


        /* ----------- Find the leading edge 50% mass value ----------- */
        while (currPtr >= &inMSDataList->mass[0]) 
        {
            if (currPtr->intensity < halfIntensity) break;
            currPtr--;
        }
        
        /* Use the two points that flank the half-intensity value to
           calculate the mass at exactly half intensity. */
        slope = (REAL_4)(((currPtr + 1)->intensity - currPtr->intensity)/
                                                 ((currPtr + 1)->mOverZ - currPtr->mOverZ));
                                                 
        intercept = (REAL_4)((currPtr->intensity) - (currPtr->mOverZ * slope));
        
        leadingMass = (REAL_4)((halfIntensity - intercept)/slope);
        
        /* ----------- Find the trailing edge 50% mass value ----------- */
        currPtr = &inMSDataList->mass[topIndex];

        while (currPtr < ptrOfNoReturn) 
        {
            if (currPtr->intensity < halfIntensity) break;
            currPtr++;
        }
        
        /* Use the two points that flank the half-intensity value to
           calculate the mass at exactly half intensity. */
        slope = (REAL_4)((currPtr->intensity - (currPtr - 1)->intensity)/
                                                 (currPtr->mOverZ - (currPtr - 1)->mOverZ));
                                                 
        intercept = (REAL_4)(((currPtr - 1)->intensity) - ((currPtr - 1)->mOverZ * slope));
        
        trailingMass = (REAL_4)((halfIntensity - intercept)/slope);
        
        /* XXXXXXXXX Why multiply by two? - JAT */
        peakWidth = (REAL_4)(trailingMass - leadingMass) * 2;
        peakWidthSum += peakWidth;
        peakWidthSquaredSum += peakWidth * peakWidth;

        /* Replace the mass with the peak width (for throwing out outliers) */
        bigTreeList->mass[i].mOverZ = (REAL_4)(trailingMass - leadingMass) * 2;        
    }


    if (bigTreeList->numObjects > 0)
    {
        /* Calculate the average peak width for the big trees. */
        avgPeakWidth = (REAL_4)(peakWidthSum/bigTreeList->numObjects);

        /* Calculate the standard deviation. */
        stdDev = sqrt((peakWidthSquaredSum/bigTreeList->numObjects) - 
                      (avgPeakWidth * avgPeakWidth));
                      
        /* Throw away peaks too far away from the average. */
        tolerance = 1.5 * stdDev;
        peakWidthSum = 0;
        for (i = 0; i < bigTreeList->numObjects; i++)
        {    
            /* Remember that the big tree mass in now really the peak width */
            if (bigTreeList->mass[i].mOverZ < (avgPeakWidth - tolerance) 
                || bigTreeList->mass[i].mOverZ > (avgPeakWidth + tolerance))
            {
                RemoveFromList(i, bigTreeList);
                i--;
            }
            else
            {
                peakWidthSum += bigTreeList->mass[i].mOverZ;
            }
        }
        avgPeakWidth = (REAL_4)(peakWidthSum/bigTreeList->numObjects);
    }
    else 
    {
        avgPeakWidth = 3;
    }
    
    if(avgPeakWidth >= 2.5)    
    {
        /*For the broad peaks, I force the fragment error to be at least 1 Da.*/
        if(gParam.fragmentErr < 1.0)
        {
            gParam.fragmentErr = 1;
        }
    }
    else
    {
        if(avgPeakWidth >= 1.5)
        {
            if(gParam.fragmentErr < 0.75)
            {
                gParam.fragmentErr = 0.75;
            }
        }
        else
        {
            /* XXXXXXX Should this be done when data is high res? JAT */
            if(gParam.fragmentErr < 0.5)
            {
                gParam.fragmentErr = 0.5;
            }
        }
    }
    
    
    avgPeakWidth = avgPeakWidth / 2;    /* the program anticipates a number that is half of the 
                                           actual peak width */

    
    if (bigTreeList) DisposeList(bigTreeList);

    return(avgPeakWidth);

}

/*
//--------------------------------------------------------------------------------
//  IntensityDescendSortFunc()
//--------------------------------------------------------------------------------
//  Modified 03.13.00 JAT - Added the secondary mass key to eliminate problems with
//                          platform specific differences.
//  
*/
INT_4 IntensityDescendSortFunc(const void *n1, const void *n2) 
{

    tMSData *n3, *n4;
    
    n3 = (tMSData *)n1;
    n4 = (tMSData *)n2;
    
    if (n3->intensity != n4->intensity) 
    {
        return (INT_4)(n3->intensity < n4->intensity)? 1:-1;
    }
    else 
    {
        if (n3->mOverZ != n4->mOverZ) 
        {
            return (INT_4)(n3->mOverZ < n4->mOverZ)? 1:-1;
        }
        return 0;
    }        
}

/*
//--------------------------------------------------------------------------------
//  MassAscendSortFunc()
//--------------------------------------------------------------------------------
*/
INT_4 MassAscendSortFunc(const void *n1, const void *n2) 
{

    tMSData *n3, *n4;
    
    n3 = (tMSData *)n1;
    n4 = (tMSData *)n2;
    
    if(n3->mOverZ != n4->mOverZ) 
    {
        return (INT_4)(n3->mOverZ > n4->mOverZ)? 1:-1;
    }
    else 
    {
        return 0;
    }        
}


/***********************************RemoveIsotopes***************************************
*
*    This function removes peaks that differ by one dalton and appear to be due to 
*   the presence of isotopes.  The algorithm is very simple and basic and only 
*   worries about whether the next ion up is an isotopic peak; it doesn't worry 
*   about two daltons up cuz I'm assuming that those ions will weeded out based on 
*   intensity.  This function won't be suitable for use on high resolution data 
*   obtained from, say, a Q-TOF.  I'll burn that bridge when I get to it.
*/

void RemoveIsotopes(tMSDataList *inMSDataList)
{
    tMSData     *currPtr = NULL;
    tMSData     *isotopePtr = NULL;
    tMSData     *isotope2Ptr = NULL;
    tMSData     *ptrOfNoReturn = NULL;
    REAL_4         upperLimit;
    REAL_4        lowerLimit;
    REAL_4        massDiff;
    REAL_4        obsdIntensityRatio;
    REAL_4        calcIntensityRatio;
    
    /* Use a tighter error since the difference is relative. */
    upperLimit = 1 + (gParam.fragmentErr / 2);
    lowerLimit = 1 - (gParam.fragmentErr / 2);
    
    currPtr     = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    while (currPtr < ptrOfNoReturn - 1) 
    {
        isotopePtr = currPtr + 1;
        while (isotopePtr < ptrOfNoReturn
               && isotopePtr->mOverZ <= currPtr->mOverZ + upperLimit) 
        {
            massDiff = isotopePtr->mOverZ - currPtr->mOverZ;
            /* Are the peaks 1 Da apart? */
            if (massDiff < upperLimit && massDiff > lowerLimit)    
            {
                /* Calculate the theoretical isotope ratio.
                *  (The 0.2 is a fudge factor.)
                */
                calcIntensityRatio = ((currPtr->mOverZ) / 1800) + 0.2;
                
                obsdIntensityRatio = (REAL_4)(isotopePtr->intensity) / (REAL_4)(currPtr->intensity);
                
                /* Does a comparison of the intensities make it look like an isotope peak? 
                *  Give 25% leeway.
                */
                if (obsdIntensityRatio <= calcIntensityRatio)/* + 0.25
                    && obsdIntensityRatio >= calcIntensityRatio - 0.25)     Fixed by RSJ*/
                {
                    /* We found what looks a +1 isotope peak, is there a peak that
                    *  looks like a +2 isotope?
                    */
                    isotope2Ptr = isotopePtr + 1;
                    while (isotope2Ptr < ptrOfNoReturn
                              && isotope2Ptr->mOverZ <= isotopePtr->mOverZ + upperLimit) 
                    {
                        massDiff = isotope2Ptr->mOverZ - isotopePtr->mOverZ;
                        /* Are the peaks 1 Da apart? */
                        if (massDiff < upperLimit && massDiff > lowerLimit)    
                        {
                            /* XXXXXX What should the equation for this ratio be? - JAT */
                            /* Calculate the theoretical isotope ratio.
                            *  (The 0.2 is a fudge factor.)
                            */
                            calcIntensityRatio = ((isotopePtr->mOverZ) / 1800) + 0.2;
                            
                            obsdIntensityRatio = (REAL_4)(isotope2Ptr->intensity) / (REAL_4)(isotopePtr->intensity);
                            
                            /* Does a comparison of the intensities make it look like an isotope peak? 
                            *  Give 25% leeway.
                            */
                            if (obsdIntensityRatio <= calcIntensityRatio)/* + 0.25
                    && obsdIntensityRatio >= calcIntensityRatio - 0.25)     Fixed by RSJ*/    
                            {
                                /* We found what looks a +2 isotope peak, Whack it. */
                                RemoveFromList((isotope2Ptr - &inMSDataList->mass[0]), inMSDataList);
                                ptrOfNoReturn--;
                                isotope2Ptr--;
                            }
                        }
                        isotope2Ptr++;                
                    }
                    
                    /* Whack the +1 isotope. */
                    RemoveFromList((isotopePtr - &inMSDataList->mass[0]), inMSDataList);
                    ptrOfNoReturn--;
                    isotopePtr--;
                }
            }
            isotopePtr++;
        }    
    
        currPtr++;
    }    
    
    return;
}

/***********************************FindMedian**********************************************
*
*    FindMedian finds the median threshold value.
*/

INT_4 FindMedian(struct MSData *firstPtr)
{
    
    struct MSData *currPtr, *biggestPtr;
    
    INT_4 biggestIntensity, lowIntensityValue, highIntensityValue;
    INT_4 numberOfIons, targetNumOfIons, median, signal;
        
/*    Count the number of datapoints, and divide in half.*/
    targetNumOfIons = 0;
    currPtr = firstPtr;
    while(currPtr != NULL)
    {
        targetNumOfIons += 1;
        currPtr = currPtr->next;
    }
    targetNumOfIons = targetNumOfIons / 2;
    
/*    Make the highest intensity datapoints negative.*/
    numberOfIons = 0;
    while(numberOfIons <= targetNumOfIons)
    {
        currPtr = firstPtr;
        biggestPtr = currPtr;
        biggestIntensity = currPtr->intensity;
        numberOfIons += 1;
        while(currPtr != NULL)
        {
            if(currPtr->intensity > biggestIntensity)
            {
                biggestIntensity = currPtr->intensity;
                biggestPtr = currPtr;
            }
            currPtr = currPtr->next;
        }
        biggestPtr->intensity = biggestPtr->intensity * -1;
    }
    
/*    Find the highest intensity datapoint from the low intensity half of the set.*/
    currPtr = firstPtr;
    biggestIntensity = currPtr->intensity;
    while(currPtr != NULL)
    {
        if(currPtr->intensity > biggestIntensity)
        {
            biggestIntensity = currPtr->intensity;
        }
        currPtr = currPtr->next;
    }
    lowIntensityValue = biggestIntensity;

/*    Find the highest intensity datapoint from the high intensity half of the set.  Recall
    that the high intensity half of the data set is negative, so I'll actually be finding
    the lowest intensity datapoint.*/
    currPtr = firstPtr;
    while(currPtr != NULL)    /*First find a negative intensity.*/
    {
        if(currPtr->intensity < 0)
        {
            biggestIntensity = currPtr->intensity;
            break;
        }
        currPtr = currPtr->next;
    }
    
    currPtr = firstPtr;        /*Now go look for the correct value.*/
    while(currPtr != NULL)    
    {
        if(currPtr->intensity > biggestIntensity && currPtr->intensity < 0)
        {
            biggestIntensity = currPtr->intensity;
        }
        currPtr = currPtr->next;
    }
    highIntensityValue = biggestIntensity * -1;

/*    Take an average of the two datapoints, which will be the median.*/
    median = (lowIntensityValue + highIntensityValue) / 2;

/*    The signal threshold is determined from the median and the user input "ionThreshold".*/
    signal = median * gParam.ionThreshold;
    
    return(signal);
}

/************************ZeroTheIons***************************************************
*
*    This function inputs the linked list of weight averaged ions (a struct of type MSData),
*    plus the REAL_4 'peakWidth', which is one-half of the width of a peak near its base.
*    It finds ions that are too close together (less than 'peakWidth') and zero's the intensity
*    field of the ion with the lowest intensity.  Those ions w/ zero intensity are free'd and
*    the list is re-linked.  It returns a pointer to this list of structs of type MSData.
*/

struct MSData *ZeroTheIons(struct MSData *firstAvMassPtr)
{
    struct MSData *currPtr, *nextPtr, *previousPtr, *structToFreePtr;
    REAL_8 diff;

    currPtr = firstAvMassPtr;    
    
    while(currPtr != NULL)
    {
        nextPtr = currPtr->next;
        
        while(nextPtr != NULL)
        {
            diff = fabs((nextPtr->mOverZ) - (currPtr->mOverZ));

            if(diff <= gParam.peakWidth * 1.5)    /* *1.75 was empirically derived*/
            {
                if(currPtr->intensity < nextPtr->intensity)    /*Which ion should be zeroed?*/
                {
                    currPtr->intensity = 0;
                }
                else
                {
                    nextPtr->intensity = 0;
                }
            }
            nextPtr = nextPtr->next;
        }
        
        /*    
        *    currPtr is moved up to the next position, and if the intensity of that value 
        *    has not been zeroed, then it breaks out of the while loop and becomes the next 
        *    currPtr.  If currPtr reaches the NULL value, the while loop terminates, and the 
        *    NULL currPtr also terminates the next loop up in the hierarchy.
        */
        while(currPtr != NULL)    
        {
            currPtr = currPtr->next;
            if(currPtr == NULL || currPtr->intensity != 0)
            {
                break;    /*Break out - you've found an ion value that is not NULL and has 
                        positive intensity.*/
            }
        }
    }

    /*
    *    Next I weed out the zero intensity ions, and re-link the non-zero ions.
    *    The original linked list of mass spectral data is free'ed.
    */

    currPtr = firstAvMassPtr;
    
    while(currPtr->intensity == 0)    /*Find the first ion that has a non-zero intensity.*/
    {
        structToFreePtr = currPtr;
        currPtr = currPtr->next;
        free(structToFreePtr);
    }
    
    firstAvMassPtr = currPtr;    /*Set the new firstAvMassPtr; this is the return value.*/
    previousPtr = currPtr;
    currPtr = currPtr->next;

    while(currPtr != NULL && previousPtr != NULL)
    {
        if(currPtr->intensity == 0)
        {
            previousPtr->next = currPtr->next;
            free(currPtr);
            currPtr = previousPtr->next;
        }
        else
        {
            previousPtr = currPtr;
            currPtr = currPtr->next;
        }
    }
    
    return(firstAvMassPtr);
}

/*******************************WeedTheIons****************************************
*
*    This function is called when the actual number of ions in a linked list of structs
*    of type MSData exceeds the value "finalIonCount".  The most intense ions are saved,
*    and the linked list is modified to remove the low intensity ions.  The discarded
*    structs are free'ed.
*/
void WeedTheIons(tMSDataList *inMSDataList, INT_4 finalIonCount, BOOLEAN spareGoldenBoys)
{
    tMSData     *currPtr;
    tMSData     *ptrOfNoReturn;
    BOOLEAN     thumbsDown;
    REAL_4         immonium[15];
    REAL_4        precursor;
    INT_4        i;
    
/*initialize immonium ions*/
    immonium[0] = gMonoMass[D] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[1] = gMonoMass[N] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[2] = gMonoMass[E] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[3] = gMonoMass[Q] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[4] = gMonoMass[H] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[5] = gMonoMass[L] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[6] = gMonoMass[M] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[7] = gMonoMass[F] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[8] = gMonoMass[P] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[9] = gMonoMass[S] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[10] = gMonoMass[T] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[11] = gMonoMass[W] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[12] = gMonoMass[Y] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[13] = gMonoMass[V] - gElementMass[1] - gElementMass[3] + gElementMass[0];
    immonium[14] = gMonoMass[K] - gElementMass[1] - gElementMass[3] - 2 * gElementMass[0] - gElementMass[2];

    precursor = (gParam.peptideMW + gParam.chargeState * gElementMass[HYDROGEN]) / gParam.chargeState;

    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];

    /* Increase the intensity of ions larger than the precursor as befitting
       their importance so they are less likely to be purged. */
    while (currPtr < ptrOfNoReturn) 
    {
        if (currPtr->mOverZ > precursor) 
        {
            currPtr->intensity *= 2.5;
        }
        currPtr++;
    }

    /* Sort the MSData in order of decreasing intensity */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),IntensityDescendSortFunc);

    /* Now that the ions are sorted, take the intensity increase given to
       the ions above the precursor back out. */
    currPtr = &inMSDataList->mass[0];
    while (currPtr < ptrOfNoReturn) 
    {
        if (currPtr->mOverZ > precursor) 
        {
            currPtr->intensity /= 2.5;
        }
        currPtr++;
    }

    /* Keep the flowers and remove the weeds. */
    currPtr = &inMSDataList->mass[finalIonCount];
    while (currPtr < ptrOfNoReturn) 
    {
        thumbsDown = TRUE;
        if (TRUE == spareGoldenBoys 
            && currPtr->normIntensity == 1 
            && currPtr->intensity > 0) 
        {
            /* Spare the golden boys */
            thumbsDown = FALSE;
        }
        else if (currPtr->mOverZ < 160)
        {
            /* Spare potential immonium ions */
            for(i = 0; i < 15; i++)
            {
                if (currPtr->mOverZ <= immonium[i] + gParam.fragmentErr 
                    && currPtr->mOverZ >= immonium[i] - gParam.fragmentErr)
                {
                    thumbsDown = FALSE;
                    break;
                }
            }
        }
        
        if (TRUE == thumbsDown) {
            RemoveFromList(currPtr - &inMSDataList->mass[0], inMSDataList);
            ptrOfNoReturn--;
            currPtr--;
        }
        
        currPtr++;
    }


    /* Resort the MSData in order of increasing mass */
    qsort(inMSDataList->mass,(size_t)inMSDataList->numObjects,
          (size_t)sizeof(tMSData),MassAscendSortFunc);


    /* Update the index values */
    for (i = 0; i < inMSDataList->numObjects; i++)
    {
        inMSDataList->mass[i].index = i;
    }


    return;
}
/****************************countIons*********************************************
*    This function counts the number of ions in the linked list of structs of type 
*    MSData.  It returns a INT_4 corresponding to the number of ions counted.
*/
INT_4 countIons(struct MSData *firstAvMassPtr)
{
    struct MSData *currPtr;
    INT_4 count = 0;
    
    currPtr = firstAvMassPtr;
    while(currPtr != NULL)
    {
        if(currPtr->normIntensity == 0/* && currPtr->mOverZ > 146.5*/)
        {
            count++;
        }
        currPtr = currPtr->next;
    }
    
    return(count);
}


/****************************countIonsAgain*********************************************
*    This function counts the number of ions in the linked list of structs of type 
*    MSData.  It returns a INT_4 corresponding to the number of ions counted.
*/
INT_4 countIonsAgain(struct MSData *firstAvMassPtr)
{
    struct MSData *currPtr;
    INT_4 count = 0;
    
    currPtr = firstAvMassPtr;
    while(currPtr != NULL)
    {
        currPtr = currPtr->next;
        count++;
    }
    
    return(count);
}
    

/****************************RemovePrecursors**************************************
*
*    I don't see why I keep these ions around.  Lets get rid of them here.
*
*/
void RemovePrecursors(tMSDataList *inMSDataList)
{
    tMSData     *currPtr;
    tMSData     *ptrOfNoReturn;
    REAL_4 precursor, precurMin2W, precurMinWA, precurMin2A, precurMinW, precurMinA;
    REAL_4 tolerance;
    
    currPtr     = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];

    tolerance = gParam.fragmentErr;
    precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / 
                gParam.chargeState;
    precurMin2W = (gParam.peptideMW - WATER - WATER + 
                  (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
    precurMinWA = (gParam.peptideMW - WATER - AMMONIA +
                  (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
    precurMin2A = (gParam.peptideMW - AMMONIA - AMMONIA + 
                  (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
    precurMinW = (gParam.peptideMW - WATER + 
                 (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;
    precurMinA = (gParam.peptideMW - AMMONIA + 
                 (gParam.chargeState * gElementMass[HYDROGEN])) / gParam.chargeState;

    while (currPtr < ptrOfNoReturn) 
    {
        if((currPtr->mOverZ <= precursor + tolerance &&  currPtr->mOverZ >= precursor - tolerance) ||
            (currPtr->mOverZ <= precurMin2W + tolerance &&  currPtr->mOverZ >= precurMin2W - tolerance) ||
            (currPtr->mOverZ <= precurMinWA + tolerance &&  currPtr->mOverZ >= precurMinWA - tolerance) ||
            (currPtr->mOverZ <= precurMin2A + tolerance &&  currPtr->mOverZ >= precurMin2A - tolerance) ||
            (currPtr->mOverZ <= precurMinW + tolerance &&  currPtr->mOverZ >= precurMinW - tolerance) ||
            (currPtr->mOverZ <= precurMinA + tolerance &&  currPtr->mOverZ >= precurMinA - tolerance))
        {
            if(currPtr->intensity > 1)
            {
                currPtr->intensity = 1;    /*if there are too many ions, then the low intensity
                                        will cause this to be removed*/
            }
            else
            {
                currPtr->intensity = 0;    /*usually intensity is a large integer, but in case it
                                        ever becomes a REAL_4 between 0 and 1, I'll stick this 
                                        condition in here - it probably won't ever be used.*/
            }
        }
        currPtr++;
    }

    return;
}

/****************************EliminateBadHighMassIons******************************
*
*    This function removes high mass ions that could not be either b or y ions, via the
*    loss of one or two amino acids.  Since ions may be due to the loss of three amino acids
*    this function stops eliminating ions that could be due to the loss of one glycine
*    and two alanines.  The combination of three glycines matches one Asn and one Gly. 
*
*/

void EliminateBadHighMassIons(tMSDataList *inPeakList)
{
    tMSData     *currPtr = NULL;
    tMSData     *windowStartPtr = NULL;
    tMSData     *ptrOfNoReturn = NULL;
    REAL_4        stopSearch, bIon, yIon;
    INT_4        i, j;
    char        test;
    
    if(gParam.peptideMW > gParam.monoToAv)
        return;    /*this function designed for monoisotopic masses only*/
        
    stopSearch = gParam.peptideMW + gElementMass[HYDROGEN] - 2 * gMonoMass[A] - 
                    gMonoMass[G] + gParam.fragmentErr;
    
    /* make all indexes positive; bad ions are later given negative index values*/
    currPtr = &inPeakList->mass[0];    
    ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];    
    while(currPtr < ptrOfNoReturn)
    {
        currPtr->index = 1;
        currPtr++;
    }

    /*look at ions from high mass to low mass*/
    currPtr = &inPeakList->mass[ (inPeakList->numObjects) - 1 ];
    ptrOfNoReturn = &inPeakList->mass[0];
    
    while(currPtr >= ptrOfNoReturn && currPtr->mOverZ > stopSearch)
    {
        test = TRUE;
        for(i = 0; i < gAminoAcidNumber; i++)
        {
            bIon = gParam.peptideMW - gParam.modifiedCTerm - gMonoMass[i];
            if(currPtr->mOverZ <= bIon + gParam.fragmentErr &&
                currPtr->mOverZ >= bIon - gParam.fragmentErr)
            {
                test = FALSE;
            }
            yIon = gParam.peptideMW + gElementMass[HYDROGEN] - gMonoMass[i];
            if(currPtr->mOverZ <= yIon + gParam.fragmentErr &&
                currPtr->mOverZ >= yIon - gParam.fragmentErr)
            {
                test = FALSE;
            }        
            for(j = 0; j < gAminoAcidNumber; j++)
            {
                bIon = gParam.peptideMW - gParam.modifiedCTerm 
                        - gMonoMass[i] - gMonoMass[j];
                            if(currPtr->mOverZ <= bIon + gParam.fragmentErr &&
                currPtr->mOverZ >= bIon - gParam.fragmentErr)
                {
                    test = FALSE;
                }
                yIon = gParam.peptideMW + gElementMass[HYDROGEN] - gMonoMass[i] - gMonoMass[j];
                if(currPtr->mOverZ <= yIon + gParam.fragmentErr &&
                currPtr->mOverZ >= yIon - gParam.fragmentErr)
                {
                    test = FALSE;
                }
            }
        }
        if(test)
        currPtr->index = -1;
        currPtr--;
    }
    
    /*get rid of the peaks with neg indexes*/
    currPtr = &inPeakList->mass[0];    
    ptrOfNoReturn = &inPeakList->mass[inPeakList->numObjects];
    while(currPtr < ptrOfNoReturn)
    {
        if(currPtr->index == -1)
        {
            RemoveFromList(currPtr - &inPeakList->mass[0], inPeakList);
            ptrOfNoReturn--;
            currPtr--;
        }

        currPtr++;
    }    
    return;
}

/****************************WindowFilter******************************************
*
*    Next the program checks to see if there are too many ions clustered together.  It does
*    this by counting the number of ions within windows of width 120 Da and making sure that
*    only a certain number of ions (ionsPerWindow) are present within any given window.  If
*    there are too many ions, it throws out those with the lowest intensity.
*/

void WindowFilter(tMSDataList *inMSDataList)
{
    tMSData     *currPtr = NULL;
    tMSData     *windowStartPtr = NULL;
    tMSData     *ptrOfNoReturn = NULL;
    INT_4         ionsInWindow;
    INT_4        endingMass;
    INT_4        charge;
    INT_4        ionsRemoved;
    REAL_4         nextChargeWindowStart;
    
        
        if(gParam.maxent3)
        {
                charge = 1;
        }
        else
        {
                charge = gParam.chargeState;    /*This is decremented to 1.*/
        }
    nextChargeWindowStart = (gParam.peptideMW + charge) / charge;
    
    
    /* Find the start of the first window at a mass greater than 176 Da. 
    *  Below this mass there is no filtering of ions.  y1 for arg is 175 
    */
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[ inMSDataList->numObjects ];
    
    while (currPtr < ptrOfNoReturn && currPtr->mOverZ < 176) currPtr++; 
    if (currPtr == ptrOfNoReturn) return;

    windowStartPtr = currPtr;    
    
    /* If the mass of the window exceeds the nextWindowStart, then the 
    *  value of charge is decremented and a nextWindowStart is calculated.  
    *  Charge can never be less than one. The variable 'nextWindowStart' 
    *  contains the m/z of the point where the window width should be changed. 
    *  For example, if the chargeState is 2 then the region below the precursor 
    *  ion has a width of SPECTRAL_WINDOW_WIDTH / 2, and the region above has 
    *  a window width of SPECTRAL_WINDOW_WIDTH / 1.
    */
    while (charge > 0 && currPtr < ptrOfNoReturn && windowStartPtr < ptrOfNoReturn) 
    {
        ionsInWindow = 0;
        
        if (currPtr->mOverZ > nextChargeWindowStart) 
        {
            charge--;
            if (charge == 0) break;
            nextChargeWindowStart = (gParam.peptideMW + charge) / charge;
        }

        endingMass = currPtr->mOverZ + (REAL_4)(SPECTRAL_WINDOW_WIDTH / charge);
    
    
        /* Ions are counted up to the endingMass value. */
        while (currPtr < ptrOfNoReturn 
               && currPtr->mOverZ < endingMass) 
        { 
            /* Only count ions that are not goldenBoys */
            if (currPtr->normIntensity == 0) ionsInWindow++;    
        
            currPtr++;
        }
        
        if (ionsInWindow > gParam.ionsPerWindow) {
            /* If there are too many ions, then purge them. */
            ionsRemoved = PurgeTheWindow(inMSDataList, windowStartPtr, ionsInWindow, endingMass);
            /*currPtr -= ionsRemoved;  */ 
            ptrOfNoReturn -= ionsRemoved;   
            /*windowStartPtr += gParam.ionsPerWindow;*/   /*Fixed by RSJ*/   
                             
        }    /*Fixed by RSJ*/ 
        /*else    /*Fixed by RSJ*/ 
        /*{    /*Fixed by RSJ*/ 
            /*windowStartPtr = currPtr;    /*Fixed by RSJ*/ 
        /*}*/
        windowStartPtr++;
        currPtr = windowStartPtr; 

    }
            
    return;
}

/****************************PurgeTheWindow****************************************
*
*    This function finds the lowest intensity ions within a particular m/z window,
*    and purges from the linked list of mass spectral CID data those ions of lowest
*    intensity.  It relinks the list and free's the space that is no longer used.
*/

INT_4 PurgeTheWindow(tMSDataList *inMSDataList, tMSData *windowStartPtr, 
                    INT_4 ionsInWindow, INT_4 endingMass)
{
    tMSData     *currPtr = NULL;
    tMSData     *ptrOfNoReturn = NULL;
    INT_4 excessIonNum, smallestIntensity, indexToRemove;
    REAL_4     precursor, precurMinW, precurMinA, precurMin2W, precurMin2A, precurMinWA;
    REAL_4  tolerance;

    
    precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN])) / 
                gParam.chargeState;
    precurMinW = (precursor - WATER) / gParam.chargeState;
    precurMinA = (precursor - AMMONIA) / gParam.chargeState;
    precurMin2W = (precursor - (2 * WATER)) / gParam.chargeState;
    precurMin2A = (precursor - (2 * AMMONIA)) / gParam.chargeState;
    precurMinWA = (precursor - WATER - AMMONIA) / gParam.chargeState;
    
    excessIonNum = ionsInWindow - gParam.ionsPerWindow;

    tolerance = gParam.fragmentErr;
    
    currPtr = windowStartPtr;
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];

    /* Whack any precursor ions first */
    if ((currPtr->mOverZ < precurMin2W - tolerance 
         && endingMass > precurMin2W - tolerance)
        || (currPtr->mOverZ < precursor + tolerance 
            && endingMass > precursor + tolerance)) 
    {

        while (currPtr < ptrOfNoReturn
               && (currPtr->mOverZ < endingMass)) 
        { 
            if (currPtr->mOverZ > precurMin2W - tolerance) {
                if ((currPtr->mOverZ >= precursor - tolerance   && currPtr->mOverZ <= precursor + tolerance)
                    || (currPtr->mOverZ >= precurMinW - tolerance  && currPtr->mOverZ <= precurMinW + tolerance)
                    || (currPtr->mOverZ >= precurMinA - tolerance  && currPtr->mOverZ <= precurMinA + tolerance)
                    || (currPtr->mOverZ >= precurMin2W - tolerance && currPtr->mOverZ <= precurMin2W + tolerance)
                    || (currPtr->mOverZ >= precurMin2A - tolerance && currPtr->mOverZ <= precurMin2A + tolerance)
                    || (currPtr->mOverZ >= precurMinWA - tolerance && currPtr->mOverZ <= precurMinWA + tolerance))
                {
                    RemoveFromList(currPtr - &inMSDataList->mass[0], inMSDataList);
                    ptrOfNoReturn--;
                    currPtr--;
                    excessIonNum--;
                }
            }
            currPtr++;
        }    
    }

    /* Now shed remaining excess ions. */
    while (excessIonNum > 0) { 
        currPtr = windowStartPtr;
        ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];

        /* Don't wipe a golden boy */
        while (currPtr->normIntensity == 1
               && currPtr < ptrOfNoReturn) currPtr++;
        if (currPtr == ptrOfNoReturn 
            || currPtr->mOverZ > endingMass) break;
               
        smallestIntensity = currPtr->intensity;
        indexToRemove = currPtr - &inMSDataList->mass[0];
        
        while (currPtr < ptrOfNoReturn
               && (currPtr->mOverZ < endingMass)) 
        { 
            if (currPtr->normIntensity != 1
                && currPtr->intensity < smallestIntensity)
            {
                smallestIntensity = currPtr->intensity;
                indexToRemove = currPtr - &inMSDataList->mass[0];
            }
            currPtr++;
        }
        RemoveFromList(indexToRemove, inMSDataList);
        excessIonNum--;
    }

    /* Return the number of ions removed. */
    return((ionsInWindow - gParam.ionsPerWindow) - excessIonNum);


}    
    
    
/***********************AddTheIonOffset********************************************
*
*    This function adds ionOffset (a user controlled variable) to the m/z field of
*    a list of structs of type MSData.
*
*/

void AddTheIonOffset(tMSDataList *inMSDataList)
{
    tMSData         *currPtr;
    tMSData         *ptrOfNoReturn;
    
    
    if (gParam.ionOffset == 0) return;
    
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    while (currPtr < ptrOfNoReturn) 
    {
        currPtr->mOverZ = currPtr->mOverZ + gParam.ionOffset;
        currPtr++;
    }

    return;
}


/***********************CheckTheIntensity******************************************
*
*    This function starts adding up the total intensity of the ions in the linked
*    list (which starts with the pointer to a struct of type MSData called firstAvMassPtr).
*    If the sum exceeds 2 billion, then all of the ions are attenuated ten fold.
*/
void CheckTheIntensity(tMSDataList *inMSDataList)
{
    INT_4             intensitySum = 0;
    tMSData         *currPtr;
    tMSData         *ptrOfNoReturn;


    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    while (currPtr < ptrOfNoReturn) 
    {
        intensitySum += currPtr->intensity;
        if(intensitySum > 2000000000)
        {
            currPtr = &inMSDataList->mass[0];
            while (currPtr < ptrOfNoReturn) 
            {
                currPtr->intensity = currPtr->intensity / 10;
                currPtr++;
            }
            currPtr = &inMSDataList->mass[0];    /*Reset the currPtr to the beginning.*/
            intensitySum = 0;    /*Reinitialize the summed intensity to zero.*/
            continue;    /*Break out of the current loop, but continue using the while loop.*/
        }
        currPtr++;    /*If nothing happens, then move on to the next ion.*/
    }

    return;
}

/***********************SortByMass*************************************************
*
*    This function takes in a pointer to the first element in the linked list of
*    structures of type MSData that contains the smoothed and average mass values
*    plus intensities.  It sorts the contents by mass - the lowest mass ions
*    are placed first in the list.
*/

void SortByMass(struct MSData *firstAvMassPtr)
{
    struct MSData *currPtr, *testPtr;
    INT_4 lowMassIntensity, highMassIntensity;
    char test;
    REAL_4 lowMOverZ, highMOverZ;
    
    currPtr = firstAvMassPtr;    /*Two pointers are compared, currPtr is supposed to contain
                                    the lowest mass ion.*/
    
    while(currPtr != NULL)    /*If currPtr is NULL then the end of the list has been reached.*/
    {
        testPtr = currPtr->next;    /*testPtr contains the high mass ion.*/
        test = TRUE;
        while(testPtr != NULL)
        {
            if(testPtr->mOverZ < currPtr->mOverZ)    /*Make the comparison.*/
            {
                lowMassIntensity = testPtr->intensity;    /*Store ion values in these INT_4 ints.*/
                lowMOverZ = testPtr->mOverZ;
                highMassIntensity = currPtr->intensity;
                highMOverZ = currPtr->mOverZ;
                
                currPtr->intensity = lowMassIntensity;    /*Swap the values back into the pointers.*/
                currPtr->mOverZ = lowMOverZ;
                testPtr->intensity = highMassIntensity;
                testPtr->mOverZ = highMOverZ;
                
                test = FALSE;    /*Do not allow the currPtr to be incremented.*/
                break;    /*Break to the outer loop.*/
            }
            testPtr = testPtr->next;    /*Get set to test the next ion.*/
        }
        if(test)    /*If no testPtr ion was found to be lower in mass, then test = FALSE.
                        Proceed to the next ion up.  If there was a rearrangement, then I
                        need to retest the new value that was swapped.*/
        {
            currPtr = currPtr->next;
        }
    }
    
    return;
}
/***********************SmoothCID**************************************
*
*    SmoothCID is a five point digital filter that uses Finnigan's coefficients,
*    which are 13, 27, 37, 27, and 13.  The first two and last two data points are
*    not smoothed.
*
*/

void SmoothCID(tMSDataList *inMSDataList)
{
    tMSData     *currPtr;
    tMSData     *ptrOfNoReturn;

    tMSData *ptrOne, *ptrTwo, *ptrThree, *ptrFour, *ptrFive;    /*These are the five points.*/
    INT_4 smoothedDataPoint;    /*This INT_4 is used to calculated the smoothed data point.*/
    
    
    if (inMSDataList->numObjects < 6) return;
    
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
/*    Initialize the five points. */
    ptrOne   = currPtr;
    ptrTwo   = currPtr + 1;
    ptrThree = currPtr + 2;
    ptrFour  = currPtr + 3;
    ptrFive  = currPtr + 4;
    
    currPtr  = currPtr + 5;
    
/*    Do the five point average.  currPtr always points to the next data point 
*    following the five
*    that are currently being smoothed.  If currPtr is NULL then the smoothing is stopped.
*/
    
    while(currPtr < ptrOfNoReturn)
    {
        if((ptrOne->intensity > 7000000) || (ptrTwo->intensity > 7000000) ||
            (ptrThree->intensity > 7000000) || (ptrFour->intensity > 7000000) ||
            (ptrFive->intensity > 7000000))    /*if intensity exceeds INT_4 int*/
        {
            ptrOne->intensity = ptrOne->intensity * 0.001;
            ptrTwo->intensity = ptrTwo->intensity * 0.001;
            ptrThree->intensity = ptrThree->intensity * 0.001;
            ptrFour->intensity = ptrFour->intensity * 0.001;
            ptrFive->intensity = ptrFive->intensity * 0.001;
            smoothedDataPoint = (13 * ptrOne->intensity) + (27 * ptrTwo->intensity) + 
            (37 * ptrThree->intensity) + (27 * ptrFour->intensity) + (13 * ptrFive->intensity);
            smoothedDataPoint = smoothedDataPoint / 117;
            smoothedDataPoint = smoothedDataPoint * 1000;
            ptrOne->intensity = ptrOne->intensity * 1000;
            ptrTwo->intensity = ptrTwo->intensity * 1000;
            ptrThree->intensity = ptrThree->intensity * 1000;
            ptrFour->intensity = ptrFour->intensity * 1000;
            ptrFive->intensity = ptrFive->intensity * 1000;
        }
        else
        {
            smoothedDataPoint = (13 * ptrOne->intensity) + (27 * ptrTwo->intensity) + 
            (37 * ptrThree->intensity) + (27 * ptrFour->intensity) + (13 * ptrFive->intensity);
            smoothedDataPoint = smoothedDataPoint / 117;
        }
        
        ptrThree->intensity = smoothedDataPoint;
        
        ptrOne = ptrTwo;        /*    Shift the five points up by one position.*/
        ptrTwo = ptrThree;
        ptrThree = ptrFour;
        ptrFour = ptrFive;
        ptrFive = currPtr;
        currPtr++;    
    }
    
    return;

}


/**********************FreeAllMSData***************************************************
*
*    An alternate way to free the data in a linked list.  Here the structs are free'ed
*    in non-reverse order.  I'll see if this crashes.
*/

void FreeAllMSData(struct MSData *currPtr)
{
    struct MSData *freeMePtr;
    
    
    while(currPtr != NULL)
    {
        freeMePtr = currPtr;
        currPtr = currPtr->next;
        free(freeMePtr);
    }
    return;
}


/******************LoadMSDataStruct********************************************
*
* LoadStruct puts mass to charge and intensity values into the appropriate struct field
* and returns the pointer to that struct.
*
*/
struct MSData *LoadMSDataStruct(REAL_4 massValue, INT_4 ionIntensity)
{
    struct MSData *currPtr = NULL;    


    currPtr = (struct MSData *)malloc(sizeof(struct MSData));
    if(currPtr == NULL)
    {
        printf("LoadMSDataStruct:  Out of memory\n");
        exit(1);
    }
    currPtr->mOverZ = massValue;
    currPtr->intensity = ionIntensity;
    currPtr->normIntensity = 0;
    currPtr->next = NULL;

    return(currPtr);
}

/****************AddToListNoNull**********************************************************
*
*     AddToListNoNull adds m/z and intensity values in an MSData struct to a linked list.  
* The first parameter is a pointer to a struct
* of type MSData containing the firt bit of data in the list.  The second parameter 
* is a pointer to a struct of type MSData
* containing a new piece of data to be added to the list.  The end of this linked 
* list is not signaled by the presence of a zero in the next field.
*
*/

struct MSData *AddToListNoNull(struct MSData *firstPtr, struct MSData *currPtr)
{
    struct MSData *lastPtr;
    
    lastPtr = firstPtr;
    
    if(firstPtr == NULL)
        firstPtr = currPtr;
    else
    {
        while(lastPtr->next != NULL)
            lastPtr = lastPtr->next;
        lastPtr->next = currPtr;
    }
    
    return(firstPtr);
}

/****************AddToCIDList**********************************************************
*
* AddToList adds m/z and intensity values in and MSData struct to a linked list.  
* The first parameter is a pointer to a struct
* of type MSData containing the firt bit of data in the list.  The second parameter 
* is a pointer to a struct of type MSData
* containing a new piece of data to be added to the list.  The end of this linked 
* list is signaled by the presence of a zero in the 
* next field.
*
*/

struct MSData *AddToCIDList(struct MSData *firstPtr, struct MSData *currPtr)
{
    
    if(firstPtr == NULL)
        {
        firstPtr = currPtr;
    }
        else
    {
        gLastDataPtr->next = currPtr;
    }
    
    gLastDataPtr = currPtr;
    return(firstPtr);
}

/***************ModifyList**********************************************************
* Replaces most recently added element in a linked list with a new element.  The 
* first paramter is a pointer to a struct of
* type MSData containing the first bit of data in the list.  The second paramter 
* is a pointer to a struct of type MSData
* containing a new piece of data to replace the last element in the list.   NULL 
* is placed in the next field of this struct element.
*
*/

void ModifyList(struct MSData *firstPtr, struct MSData *currPtr)
{
    struct MSData *lastPtr, *nextToLastPtr;
    
    lastPtr = firstPtr;
    
    if(firstPtr == NULL)
    {
        firstPtr = currPtr;
        return;
    }
    if(firstPtr->next == NULL)
    {
        return;
    }
    else
    {
        while(lastPtr->next != NULL)
        {
            nextToLastPtr = lastPtr;
            lastPtr = lastPtr->next;
        }
        free(lastPtr);
        nextToLastPtr->next = currPtr;
    }
    
    currPtr->next = NULL;
    gLastDataPtr = currPtr;
}




/***************FindThreshold*****************************************
*
* FindThreshold sums the intensities of all of the ions in the list  
* of tMSData structs whose intensity is greater than zero.  This sum 
* is then divided by the total number of non-zero ions present, and is 
* returned as a INT_4.
*
******/
INT_4 FindThreshold(tMSDataList *inMSDataList)
{
    tMSData     *currPtr;
    tMSData     *ptrOfNoReturn;
    INT_4         signal = 0;
    INT_4         numOfNonZeroPts = 0;
    REAL_4         intensitySum = 0;
    REAL_4        noise = 0;
    
            
    if (inMSDataList->numObjects != 0) 
    {
        currPtr = &inMSDataList->mass[0];
        ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
        
        while (currPtr < ptrOfNoReturn) 
        {
            if(currPtr->intensity != 0)    /*dont count the zero intensity data points*/
            {
                numOfNonZeroPts++;
                intensitySum += currPtr->intensity;    /*Add the intensity.*/
            }
            currPtr++;    /*Increment to the next data point.*/
        }
        
        noise = intensitySum / numOfNonZeroPts;
        printf("Average signal (noise) = %.2f\n",noise);  /*debug*/
        printf("Data points = %ld\n", numOfNonZeroPts); /*debug*/
        signal = noise * gParam.ionThreshold;
        
        if(signal == 0)
            signal = 1;    /*added 11/10/98 for qtof data*/
    }
    
    return(signal);
}
/***********************************LowMassIonRemoval****************************************
*
*    Low mass ions that could not be due to amino acid immonium ions and other small fragments
*    derived from amino acids are removed.  For tryptic peptides, the mass range below y1 for
*    Lys (m/z 147) are removed, and for other cleavages, ions below m/z 90 (y1 for Ala) are
*    eliminated.
*
*/

void    LowMassIonRemoval(tMSDataList *inMSDataList)
{
    REAL_4 lowMassIons[23] =     {44.0500,  70.0657,  87.0922, 112.0875,  87.0558,
                              88.0399, 102.0555,  84.0450, 101.0715, 129.0664,
                             110.0718,  86.0970,  86.0970,  84.0814, 101.1079, 
                             129.1028, 104.0534, 120.0813,  70.0657,  60.0449,
                              74.0606, 136.0762,  72.0813
    };
    
    INT_4 i;
    tMSData     *currPtr;
    tMSData     *ptrOfNoReturn;
    BOOLEAN     thumbsDown;
    
    currPtr = &inMSDataList->mass[0];
    ptrOfNoReturn = &inMSDataList->mass[inMSDataList->numObjects];
    
    while (currPtr < ptrOfNoReturn) 
    {
        if(currPtr->mOverZ < 74.5)
        {
            thumbsDown = TRUE;
            /* Spare potential immonium ions */
            for(i = 0; i < 23; i++)
            {
                if (currPtr->mOverZ <= lowMassIons[i] + gParam.fragmentErr &&
                    currPtr->mOverZ >= lowMassIons[i] - gParam.fragmentErr)
                {
                    thumbsDown = FALSE;
                    break;
                }
            }
            
            if (TRUE == thumbsDown) {
                RemoveFromList(currPtr - &inMSDataList->mass[0], inMSDataList);
                ptrOfNoReturn--;
                currPtr--;
            }
        }
        
        currPtr++;
    }

    return;
}



/***********************************IonCondenser**********************************************
*
*    This function was originally designed with centroided data in mind, but it seemed to work
*    well for profile data, too.  The idea is that groups of ions that are above the 'threshold' 
*    of ion intensity, and are close in m/z (as defined by 'peakWidth'). The value of 'precursor'
*   is the m/z value of the precursor ion, and it is used to determine when to switch the 
*   threshold value to one-half of what was read from the file Lutefisk.params.
*     
*/

tMSDataList *IonCondenser(tMSDataList *inMSDataList) 
{
    
    INT_4         i, j;    
    REAL_4         halfWidth = gParam.peakWidth/2;
    tMSDataList *  intensityOrderedList = NULL;
    tMSDataList *  peakList = NULL;
    tMSData     potentialPeak;
    tMSData     peak;
    tMSData     msData;
    REAL_8        intensitySum = 0.0;
    INT_4         prevIntensity = 0;
    REAL_4        avgIonMass = 0;


        
    /* Add index values */    
    for (i = 0; i < inMSDataList->numObjects; i++)
    {
        inMSDataList->mass[i].index = i;
    }    

    intensityOrderedList = (tMSDataList *) CopyList( inMSDataList );
    if (!intensityOrderedList) 
    {
        printf("Ran out of memory in IonCondenser()!\n");
        exit(1);
    }
    
    /* Sort the intensityOrderedList in order of decreasing intensity */
    qsort(intensityOrderedList->mass,(size_t)intensityOrderedList->numObjects,
          (size_t)sizeof(tMSData),IntensityDescendSortFunc);


    peakList = (tMSDataList *) CreateNewList( sizeof(tMSData), 1000, 1000 );
    if (!peakList) 
    {
        printf("Ran out of memory in IonCondenser()!\n");
        exit(1);
    }


    /* Find all the peak tops that are at least a peakWidth apart. */
    for (i = 0; i < intensityOrderedList->numObjects; i++)
    {
        potentialPeak = intensityOrderedList->mass[i];
        
        /* Quit when we hit the zero intensity points */
        if (potentialPeak.intensity <= 0) break;
        
        /*    Get rid of ions below mass 43 (44 is the immonium ion of alanine). */
        if (potentialPeak.mOverZ < 43.0) continue;

        /* Don't try to make a peak out of an ion whose intensity 
           is less than the threshold. */           
        if (potentialPeak.intensity < gParam.intThreshold) 
        {
            if((potentialPeak.mOverZ > 43.0 && potentialPeak.mOverZ < 147.5) 
               || (potentialPeak.mOverZ > 174.5 && potentialPeak.mOverZ < 175.5) 
               || (potentialPeak.mOverZ > 158.5 && potentialPeak.mOverZ < 159.5))
            {
                /* Let immoniums thru even if they are low intensity */
            }
            else continue;
        }

        /* Would the potential peak overlap the domain of an existing peak? */
        for (j = 0; j < peakList->numObjects; j++)
        {
            peak = peakList->mass[j];
            
            if (potentialPeak.mOverZ >= peak.mOverZ - gParam.peakWidth
                && potentialPeak.mOverZ <= peak.mOverZ + gParam.peakWidth) break;
                
            /**Fixed by RSJ.  Try to keep peaks closer to 1 Da apart.**/
            if(gParam.peakWidth < 0.8)
            {
                if (potentialPeak.mOverZ >= peak.mOverZ - 0.8
                && potentialPeak.mOverZ <= peak.mOverZ + 0.8) break;
            }
            
        }
        if (j == peakList->numObjects) 
        {
             /* Add a new peak. */
             
             intensitySum = 0;
             avgIonMass = 0;
             
             /* Calculate weighted average mass for the peak. */

            /* Leading edge... */
            prevIntensity = potentialPeak.intensity;
            j = potentialPeak.index;
            while (j >= 0) 
            {
                msData = inMSDataList->mass[j];
                if (msData.mOverZ < potentialPeak.mOverZ - halfWidth) break;
                if (msData.intensity >= gParam.intThreshold) 
                {
                    intensitySum += msData.intensity;
                    avgIonMass   += (msData.intensity * msData.mOverZ);
                }
                else if (prevIntensity >= msData.intensity
                         && ((msData.mOverZ > 43.0 && msData.mOverZ < 147.5)
                             || (msData.mOverZ > 174.5 && msData.mOverZ < 175.5)
                             || (msData.mOverZ > 158.5 && msData.mOverZ < 159.5)))
                {
/*              else if (prevIntensity >= msData.intensity
                         || (msData.mOverZ > 43.0 && msData.mOverZ < 147.5)
                         || (msData.mOverZ > 174.5 && msData.mOverZ < 175.5)
                         || (msData.mOverZ > 158.5 && msData.mOverZ < 159.5))
                {
*/                  intensitySum += msData.intensity;
                    avgIonMass   += (msData.intensity * msData.mOverZ);                
                }
                prevIntensity = msData.intensity;
                j--;
            }
            /* Trailing edge... */
            prevIntensity = potentialPeak.intensity;
            j = potentialPeak.index + 1;
            while (j < inMSDataList->numObjects) 
            {
                msData = inMSDataList->mass[j];
                if (msData.mOverZ > potentialPeak.mOverZ + halfWidth) break;
                if (msData.intensity >= gParam.intThreshold) {
                    intensitySum += (INT_4)msData.intensity;
                    avgIonMass += (REAL_4)(msData.intensity * msData.mOverZ);
                }
                else if (prevIntensity >= msData.intensity
                         && ((msData.mOverZ > 43.0 && msData.mOverZ < 147.5)
                             || (msData.mOverZ > 174.5 && msData.mOverZ < 175.5)
                             || (msData.mOverZ > 158.5 && msData.mOverZ < 159.5)))
                {
/*              else if (prevIntensity >= msData.intensity
                         || (msData.mOverZ > 43.0 && msData.mOverZ < 147.5)
                         || (msData.mOverZ > 174.5 && msData.mOverZ < 175.5)
                         || (msData.mOverZ > 158.5 && msData.mOverZ < 159.5))
                {
*/                  intensitySum += msData.intensity;
                    avgIonMass   += (msData.intensity * msData.mOverZ);                
                }
                prevIntensity = msData.intensity;
                j++;
            }
            
            if (intensitySum == 0) {
                printf("ionCondenser(): intensitySum = 0.\n" 
                       "Threshold: %d, potential peak: %6.2f  %d (index %d)\n", 
                           gParam.intThreshold, potentialPeak.mOverZ, 
                           potentialPeak.intensity, potentialPeak.index);
                           
                /* Leading edge... */
                prevIntensity = potentialPeak.intensity;
                j = potentialPeak.index;
                while (j >= 0) 
                {
                        msData = inMSDataList->mass[j];
                        printf("Leading: %6.2f  %d (%d)\n", msData.mOverZ, msData.intensity, msData.index);
                        if (msData.mOverZ < potentialPeak.mOverZ - halfWidth) break;
                        j--;
                    }      
                /* Trailing edge... */
                    prevIntensity = potentialPeak.intensity;
                    j = potentialPeak.index + 1;
                    while (j < inMSDataList->numObjects) 
                    {
                        msData = inMSDataList->mass[j];
                        printf("Trailing: %6.2f  %d (%d)\n", msData.mOverZ, msData.intensity, msData.index);
                        if (msData.mOverZ > potentialPeak.mOverZ + halfWidth) break;
                        j++;
                    }    
                    
 /**************/
    for (j = 0; j < inMSDataList->numObjects; j++)
    {
        inMSDataList->mass[i].index = i;
printf("AFTER: %6.2f  %d\n", inMSDataList->mass[j].mOverZ, inMSDataList->mass[j].intensity);
    }    
printf("pointer: %d\n", inMSDataList);
/***********/

                      
                exit(1);
            }
    
            potentialPeak.mOverZ = avgIonMass / intensitySum;

            /* Now double check to be sure that the mass still is not within
               a peakwidth of another peak. */
            
            for (j = 0; j < peakList->numObjects; j++)
            {
                peak = peakList->mass[j];
                
                if (potentialPeak.mOverZ >= peak.mOverZ - gParam.peakWidth
                    && potentialPeak.mOverZ <= peak.mOverZ + gParam.peakWidth) break;
                
            }
            if (j == peakList->numObjects) 
            {
                 /* Add a new peak. */
                if(!AddToList(&potentialPeak, peakList)) 
                {
                    printf("Ran out of memory in IonCondenser()!\n");
                    exit(1);
                }
            }
        }
    }
    
    if (intensityOrderedList) DisposeList(intensityOrderedList);
            
    /* Resort the MSData in order of increasing mass */
    qsort(peakList->mass,(size_t)peakList->numObjects,
          (size_t)sizeof(tMSData),MassAscendSortFunc);
            
    /* Add index values */    
    for (i = 0;    i <    peakList->numObjects; i++)
    {
        peakList->mass[i].index = i;
        peakList->mass[i].normIntensity = 0;
    }    

    return peakList;
}

/* -------------------------------------------------------------------------
//  my_fgets 
*/
char *my_fgets(char *s, INT_4 n, FILE *fp)
{
  register char *t = s;
  register long c;

  if (n < 1)
    return(NULL);

  while (--n) {
    if ((c = getc(fp)) < 0) {
      if (feof(fp) && t != s) break;
      return(NULL);
    }

    *t++ = c;
    if (c == '\n' || c == '\32') break;
    /* This is to handle stupid windows files that end each line with \r\n */
    else if (c == '\r') {
      if ((c = getc(fp)) < 0) {
        if (feof(fp) && t != s) break;
        return(NULL);
      }
      if (c != '\n') {
        /* Roll back the file pointer or we will miss characters. */
        fseek(fp, -1, SEEK_CUR);
      }
      break;
    }
  }

  *t = '\0';

  return(s);
}

#ifdef DEBUG
/**************************** DumpMassList **************************************
*
*    For debugging purposes only.
*
*/

void DumpMassList(struct MSData *firstPtr) 
{

    struct MSData     *currPtr = NULL;
    INT_4 ionNum;
    INT_4 intensity[3000];    
    REAL_4 mass[3000];        

    /* Create mass and integerMass arrays from linked list data */
    ionNum = 0;
    currPtr = firstPtr;
    while(currPtr != NULL)
    {
        mass[ionNum] = currPtr->mOverZ;
        intensity[ionNum] = currPtr->intensity;
        ionNum++;
        if(ionNum >= 3000)
            break;    
        currPtr = currPtr->next;
    }
}

/**************************** DumpMSData **************************************
*
*    For debugging purposes only.
*
*/

void DumpMSData(tMSDataList *inList) 
{

    INT_4 i;
    
    printf("\n\nMASS LIST:\n");
    for (i = 0; i < inList->numObjects; i++) 
    {
        if (i >= 3000) break;    
        printf("%ld  %7.3f  %ld\n", i, inList->mass[i].mOverZ, inList->mass[i].intensity);
    
    }
    

}

#endif
