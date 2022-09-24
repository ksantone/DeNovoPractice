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
		LutefiskXP is a program designed to aid in the interpretation of CID data of peptides.  
	The main assumptions are that the data is of reasonable quality, the N- and C-terminal
	modifications (if any) are known, and the precursor ion charge (and therefore the 
	peptide molecular weight) are known.  The ultimate goal here is to develop code that
	can utilize msms data in conjunction with ambiguous and incomplete Edman sequencing data,
	sequence tags, peptide derivatization, and protein or est database searches.  An older 
	version of LutefiskXP has been written in FORTRAN and runs on 68K Macs that have an fpu
	(1991, 39th ASMS Conference on Mass Spectrometry and Allied Topics, Nashville, TN, pp 1233-
	1234).  This is a different and improved algorithm partly inspired by Fernandez-de-Cossjo, 
	et al. (1995) CABIOS Vol. 11 No. 4 pp 427-434.  Combining this msms interpretation algorithm
	with Edman sequencing, database searches, and derivatization is entirely of my own design;
	J. Alex Taylor implemented the changes in the FASTA code (Bill Pearson, U. of VA) so that
	the LutefiskXP output can be read directly by the modified FASTA program.  In addition, there
	were a number of additional critical changes made to FASTA to make it more compatible with
	msms sequencing data.
	
		The trademark LutefiskXP was chosen at random, and is not meant to imply any similarity
	between this computer program and the partially base-hydrolzyed cod fish of the same name. 
*/

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#include "LutefiskPrototypes.h"
#include "LutefiskDefinitions.h"

REAL_4 *	spectrum1 = NULL;
REAL_4 * spectrum2 = NULL;
REAL_4 * tau = NULL;

extern UINT_4 SIZEOF_SPECTRA;

static void FastFourier(REAL_4 *data, UINT_4 nn, INT_4 isign);
static void twofft(REAL_4 data1[], REAL_4 data2[], REAL_4 fft1[], REAL_4 fft2[], UINT_4 n);
static void realft(REAL_4 data[], UINT_4 n, int isign);


/*************************************CalcNormalizedExptPeaks**********************************
*
*  It takes the peaks from the current MS/MS spectrum and for each charge state (1,2,3,4,5+) 
*  it divides the spectrum into ten equal sections and normalizes the peaks in each section 
*  to an intensity of 50.  The parent (and potential parent derivatives in that charge state) 
*  are excluded from this normalization.
*/
void CalcNormalizedExptPeaks(struct MSData *firstMassPtr) 
{

	INT_4      		segment;			/* The current segment (1/10th of the spectrum). */
	REAL_4 			tolerance;			/* The daughter ion error tolerance. */
	REAL_4 			offset; 			/* The mass offset for daughter ions. */
	REAL_4 			segmentSize;		/* The size (in m/z) of a segment */
	REAL_4 			maxIntensity;		/* The max peak intensity in a segment */
	REAL_4   		normalizingFactor;  /* The normalizing factor for a segment */
	struct MSData	*pPeak;				/* Pointer to the current MS/MS peak. */
	struct MSData	*pSegment;			/* Pointer to an MS/MS peak in the segment. */
	REAL_4 			precursor;			/* The m/z of the parent ion. */
	INT_4 		charge, segmentNum;


	charge = gParam.chargeState;
	precursor = (gParam.peptideMW + (charge * gElementMass[HYDROGEN])) / charge;
	tolerance = gParam.fragmentErr;
	offset    = 0;	/*this has already been incorporated into the list of ions*/
	segmentNum = ((msms.scanMassHigh - msms.scanMassLow) / AV_RESIDUE_MASS) + 1;
	if(segmentNum == 0)
		exit(1);
	segmentSize = (msms.scanMassHigh - msms.scanMassLow) / segmentNum;

	segment = 0;	/* Index to the segment */
	maxIntensity = 0;
	pPeak = firstMassPtr;
	while(pPeak != NULL && segment < segmentNum) {
		
		if (pPeak->mOverZ > msms.scanMassLow - gParam.fragmentErr + (segment * segmentSize)) {
			
			pSegment = pPeak;
			
			/* Find the peak with the max intensity in this segment */
			while (pPeak != NULL && pPeak->mOverZ < msms.scanMassLow + ((segment+1) * segmentSize)) 
			{
				if (pPeak->mOverZ > (precursor - 2*H2O - 2*tolerance) && 
					pPeak->mOverZ < (precursor + (tolerance * 2))) 
				{
					/* Count as max only if it is not the parent or a parent derivative. */
					if (!closeEnough((pPeak->mOverZ - (precursor + offset - H2O/charge)), tolerance * 2) &&
					    !closeEnough((pPeak->mOverZ - (precursor + offset - NH3/charge)), tolerance * 2) &&
					    !closeEnough((pPeak->mOverZ - (precursor + offset - 2*H2O/charge)), tolerance * 2) &&
					    !closeEnough((pPeak->mOverZ - (precursor + offset - 2*NH3/charge)), tolerance * 2) &&
					    !closeEnough((pPeak->mOverZ - precursor + offset), tolerance * 2))
					{ 	
						if(pPeak->intensity > maxIntensity)
						{
							maxIntensity = pPeak->intensity;
						}
					}
					else 
					{
						pPeak->normIntensity = -1; /* Flag it as a parent or derivative */
					}
				}
				else 
				{
					if (pPeak->intensity > maxIntensity) 
					{
						maxIntensity = pPeak->intensity;
					}
				}				
				pPeak = pPeak->next;

			}
			
			/* Normalize the peaks in this segment */
			if ((pPeak == NULL) || (pPeak->mOverZ > pSegment->mOverZ))
			{
				if (maxIntensity > 0) 
				{
					normalizingFactor = 50.0/maxIntensity;
				}
				else 
				{
					normalizingFactor = 1;
				}
		
				while (pSegment != NULL) 
				{
					if (pSegment->normIntensity == -1) 
					{
						/* The peak should have a normalized intensity of 0 if it is flagged
						*  as a parent or parent derivative
						*/
						pSegment->normIntensity = 0;
					}
					else 
					{
						pSegment->normIntensity = normalizingFactor * pSegment->intensity;
					}
					
					pSegment = pSegment->next;
					if(pPeak == NULL) 
					{
						if(pSegment == NULL) 
						{
							break;
						}			
					}
					else if(pSegment == NULL || pSegment->mOverZ >= pPeak->mOverZ)
					{
						break;
					}
				}	 
			}
		}
		segment++;
		maxIntensity = 0;
	}
}
/*************************************FillInSpectrum1**********************************
*
*  This function takes the normalized peak intensities from the current MS/MS spectrum
*  and creates a dummied up spectrum.
*/
void FillInSpectrum1(struct MSData *firstMassPtr) 
{

	struct MSData	*pPeak;				/* Pointer to the current MS/MS peak. */
	INT_4 massInt, lowEnd, highEnd, i;
	REAL_4 precursor;
	
	precursor = (gParam.peptideMW + (gParam.chargeState * gElementMass[HYDROGEN]))
					/ gParam.chargeState;
	lowEnd = (precursor - 35) * 2;	/*From lowEnd to highEnd, I attenuate the spectrum1 values*/
	highEnd = (precursor + 2) * 2;
	
	
	if (!spectrum1) {
		return;
	}
	
/*	
	For data with higher mass accuracy, the spectrum1 is made narrower.  Since I bin every 0.5
	Da, three bins would give a peak width of about 1.5 Da.
*/
	if(gParam.fragmentErr <= 0.75)
	{
		pPeak = firstMassPtr;
		while (pPeak != NULL) 
		{
			massInt = (((pPeak->mOverZ) * 2) + 0.5);
			if (massInt > 2 && massInt < (SIZEOF_SPECTRA - 2) ) 
			{
				spectrum1[massInt] = pPeak->normIntensity;
				if(gParam.qtofErr == 0 || gParam.qtofErr >= 0.25)
				{
					if (0.5 * pPeak->normIntensity > spectrum1[massInt - 1]) 
					{
						spectrum1[massInt - 1] = 0.75 * pPeak->normIntensity;
					}
					if(0.5 * pPeak->normIntensity > spectrum1[massInt + 1])
					{
						spectrum1[massInt + 1] = 0.75 * pPeak->normIntensity;
					}
				}
			}
			pPeak = pPeak->next;
		}
	}
/*
	For data with worse errors, I go for a peak that is 5 bins wide.
*/
	else
	{
		pPeak = firstMassPtr;
		while (pPeak != NULL) 
		{
			massInt = (((pPeak->mOverZ) * 2) + 0.5);
			if (massInt > 2 && massInt < (SIZEOF_SPECTRA - 2) ) 
			{
				spectrum1[massInt] = pPeak->normIntensity;
				if (0.5 * pPeak->normIntensity > spectrum1[massInt - 1]) 
				{
					spectrum1[massInt - 1] = 0.5 * pPeak->normIntensity;
				}
				if(0.25 * pPeak->normIntensity > spectrum1[massInt - 2])
				{
					spectrum1[massInt - 2] = 0.25 * pPeak->normIntensity;
				}
				if(0.5 * pPeak->normIntensity > spectrum1[massInt + 1])
				{
					spectrum1[massInt + 1] = 0.5 * pPeak->normIntensity;
				}
				if(0.25 * pPeak->normIntensity > 0.25 * spectrum1[massInt + 2])
				{
					spectrum1[massInt + 2] = 0.25 * pPeak->normIntensity;
				}
			}
			pPeak = pPeak->next;
		}

	}
	
	for(i = lowEnd; i < highEnd; i++)
	{
		spectrum1[i] = spectrum1[i] * 0.5;
	}
}
/*************************************SetupCrossCorrelation**********************************
*
*  Sets aside memory blocks for the two spectra that will be dummied up and padded, and for
*  the array to contain the cross-correlation results, tau.
*/
void SetupCrossCorrelation(void) 
{

	INT_4 i;

	/* The spectral arrays must be of the same power of 2.  So, if the spectra was not 
	*  acquired above m/z 2048 we will uses that as the array size because it will run
	*  much faster than if we use 4096. (The cross-correlation is done at unit resolution.)
	*/
	SIZEOF_SPECTRA = SIZEOF_SPECTRA_BIG;
	
	if (spectrum1) 
	{
		free(spectrum1);   /* Throw away the old data (if any exists) */
	}	

	spectrum1 = (REAL_4 *) malloc(SIZEOF_SPECTRA*sizeof(REAL_4));
	if (NULL == spectrum1)
	{
		printf("Not enough memory to allocate spectrum1");
		exit(1);
	}

	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		spectrum1[i] = 0;
	}
	
	
	if (spectrum2) {
		free(spectrum2);   /* Throw away the old data (if any exists) */
	}	

	spectrum2 = (REAL_4 *) malloc(SIZEOF_SPECTRA*sizeof(REAL_4));
	if (NULL == spectrum2)
	{
		printf("Not enough memory to allocate spectrum2");
		exit(1);
	}

	for(i = 0; i < SIZEOF_SPECTRA; i++)
	{
		spectrum2[i] = 0;
	}
	
	
	if (tau) {	/* tau will be the array for the results of the cross-correlation
				*  and it needs to be twice as large as the spectra.
				*/
		free(tau);   /* Throw away the old data (if any exists) */
	}	

	tau = (REAL_4 *) malloc(SIZEOF_SPECTRA * 2 * sizeof(REAL_4));
	if (NULL == tau)
	{
		printf("Not enough memory to allocate tau");
		exit(1);
	}

	for(i = 0; i < (SIZEOF_SPECTRA) * 2; i++)
	{
		tau[i] = 0;
	}
}
/*************************************CrossCorrelate**********************************
*
*/
void CrossCorrelate(REAL_4 *array1, REAL_4 *array2, UINT_4 n, REAL_4 *result) 
{

	UINT_4		i;				/* Loop index. */
	REAL_4 				temp;
	REAL_4				*workSpace;
	
	
	workSpace = (REAL_4 *) malloc((n * 2) * sizeof(REAL_4));
	for(i = 0; i < (n * 2); i++)
	{
		workSpace[i] = 0;
	}
	if (workSpace) {
		workSpace--;	/* This is done so the array is not treated as 0 based. */
		twofft(array1, array2, workSpace, result, n);
		for (i = 2; i <= n+2; i += 2) {
			result[i-1] = (workSpace[i-1] * (temp = result[i-1]) + workSpace[i] * result[i])/(n * 2);
			result[i] = (workSpace[i] * temp - workSpace[i-1] * result[i])/(n * 2);
		}
		result[2] = result[n+1];
		realft(result, n, -1);
		
		workSpace++;
		free(workSpace);
	}	
}	

/*************************************FastFourier**********************************
*
*/
static void FastFourier(REAL_4 *array, UINT_4 nn, INT_4 isign) 
{

	UINT_4 		n;
	UINT_4 		mmax;
	UINT_4 		m;
	UINT_4 		j;
	UINT_4		istep;
	UINT_4		i;
	REAL_4				tempr;
	REAL_4				tempi;
	REAL_8 				wtemp;
	REAL_8 				wr;
	REAL_8			 	wpr;
	REAL_8				wpi;
	REAL_8				wi;
	REAL_8				theta;
	
	n = nn << 1;
	j = 1;
	for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(array[j],array[i]);
			SWAP(array[j+1],array[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}	
	
	mmax = 2;
	while (n> mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959/mmax);
		wpr = -2.0 * pow(sin(0.5 * theta),2);
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1; m<mmax; m+=2) {
			for (i=m; i<=n; i+=istep) {
				j = i + mmax;
				tempr = wr * array[j] - wi * array[j+1];
				tempi = wr * array[j+1] + wi * array[j];
				array[j] = array[i] - tempr;
				array[j+1] = array[i+1] - tempi;
				array[i] += tempr;
				array[i+1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
	}
}

/*************************************twofft**********************************
*
*/
static void twofft(REAL_4 *array1, REAL_4 *array2, REAL_4 *fft1, REAL_4 *fft2, UINT_4 n) 
{

	UINT_4 		nn3;
	UINT_4 		nn2;
	UINT_4		jj;
	UINT_4		j;
	REAL_4 				rep;
	REAL_4				rem;
	REAL_4				aip;
	REAL_4				aim;
	
	
	nn3 = 1 + (nn2 = 2 + n + n);
	for (j=1, jj=2; j<=n; j++, jj+=2) {
		fft1[jj-1] = array1[j];
		fft1[jj] = array2[j];
	}
	FastFourier(fft1,n,1);
	fft2[1] = fft1[2];
	fft1[2] = fft2[2] = 0.0;
	for (j=3; j<=n+1; j+=2) {
		rep = 0.5 * (fft1[j] + fft1[nn2-j]);
		rem = 0.5 * (fft1[j] - fft1[nn2-j]);
		aip = 0.5 * (fft1[j+1] + fft1[nn3-j]);
		aim = 0.5 * (fft1[j+1] - fft1[nn3-j]);
		fft1[j] = rep;
		fft1[j+1] = aim;
		fft1[nn2-j] = rep;
		fft1[nn3-j] = -aim;
		fft2[j] = aip;
		fft2[j+1] = -rem;
		fft2[nn2-j] = aip;
		fft2[nn3-j] = rem;
	}
}	
/*************************************realft**********************************
*
*/
static void realft(REAL_4 *array, UINT_4 n, int isign) 
{

	UINT_4 i, i1, i2, i3, i4, np3;
	REAL_4  c1=0.5, c2, h1r, h1i, h2r, h2i;
	REAL_8 wr, wi, wpr, wpi, wtemp, theta;
	
	theta = 3.141592653589793/(REAL_8) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		FastFourier(array, n>>1, 1);
	}
	else {
		c2 = 0.5;
		theta = -theta;
	}
	
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr  = 1.0 + wpr;
	wi  = wpi;
	np3 = n + 3;
	for (i=2; i<=(n>>2); i++) {
		i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
		h1r = c1 * (array[i1] + array[i3]);
		h1i = c1 * (array[i2] - array[i4]);
		h2r = -c2 * (array[i2] + array[i4]);
		h2i =  c2 * (array[i1] - array[i3]);
		array[i1] =  h1r + wr * h2r - wi * h2i;
		array[i2] =  h1i + wr * h2i + wi * h2r;
		array[i3] =  h1r - wr * h2r + wi * h2i;
		array[i4] = -h1i + wr * h2i + wi * h2r;
		wr = (wtemp = wr) * wpr - wi * wpi + wr;
		wi = wi * wpr + wtemp * wpi + wi;
	
	}
	if (isign == 1) {
	
		array[1] = (h1r = array[1]) + array[2];
		array[2] = h1r - array[2];
	}
	else {
		array[1] = c1 * ((h1r = array[1]) + array[2]);
		array[2] = c1 * (h1r - array[2]);
		FastFourier(array, n>>1, -1);
	}		

}
