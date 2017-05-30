#include <stdio.h>
#include <fitsio.h>
#include <string.h>
#include "excpt.hh"
#include <math.h>
#include <vector>
#include <iostream>

void dexHeapSort(int* dex, const float* data, int n);

/*
 * enclosep - reads the 2D FITS output from snalys and normalizes the
 *            probability content and returns the enclosed probability.
 *
 * Rahman Amanullah, 2004-01-22
 *
 * This is actually only C, and almost no C++.
 *
 * The dexHeapSort function was stolen from Rob Knop's omlam code.
 *
 * TODO: - Add a key that specifies the state of the data file, i.e.
 *         if it only contains the enclosed probability.
 *       - Add an extension label to the mlextension to avoid problems
 *         if this extension for some reason is no longer the first one.
 * 
 */
int main(int argc,char* argv[]) {
  fitsfile *infptr, *outfptr;   /* FITS file pointers defined in fitsio.h */
  int status = 0, ii = 1;       /* status must always be initialized = 0  */
  int k = 0, l = 0;
  char comment[FLEN_COMMENT], *outf;

  int dim, maxdim = 8, totpix = 1, anynul = 0, *dex;
  int bitpix = 0, datatype = 0;
  long naxes[maxdim], firstpx[maxdim];

  int nkeys;
  char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */

  float *data, nullvalue = 0, totsum = 0, minval = 0, maxval = 0;

  //  try {
  if ( argc != 3 ) {
    cerr << "Usage: enclosep <infile> <outfile>\n";
    exit(5);
  }

  /* Add a name to the extension and open the file with the proper
     name, so we do not have to check wether or not it is an image.
     if ( !fits_open_file(&infptr, argv[1].[mlfunction], READONLY, &status) )
  */

  outf = (char *) malloc(strlen(argv[2]) + 5); strcat( outf, "!" );
  if ( !fits_open_file(&infptr, argv[1], READONLY, &status) &&
       !fits_create_file(&outfptr, strcat(outf, argv[2]) , &status) ) {

    fits_get_img_param(infptr, maxdim, &bitpix, &dim , naxes, &status);

    fits_create_img(outfptr, bitpix, dim, naxes, &status);

    /* copy all the user keywords (not the structural keywords) */
    fits_get_hdrspace(infptr, &nkeys, NULL, &status); 
    for (ii = 1; ii <= nkeys; ii++) {
      fits_read_record(infptr, ii, card, &status);
      if (fits_get_keyclass(card) > TYP_CMPRS_KEY)
	fits_write_record(outfptr, card, &status);
    }

    /* We assume that the bitpix = -32 */
    if ( bitpix != FLOAT_IMG ) {
      cerr << "We require that BITPIX = -32\n";
      exit(5);
    }

    /* Do some initialization before reading the file */
    for ( ii = 0; ii < dim; ii++ ) {
      totpix *= naxes[ii];
      firstpx[ii] = 1;
    }
    data = (float *) malloc(totpix * sizeof(float)); 

    /* Read the ML function */
    fits_read_pix(infptr, TFLOAT, firstpx, totpix, 
		  NULL, data , NULL, &status);
    if (status == END_OF_FILE)  status = 0; /* Reset after normal error */    
    fits_close_file(infptr, &status);

    /* Find the max and min values */
    minval = maxval = data[0];
    for ( k = 0; k < totpix; k++ ) {
      if ( data[k] < minval ) { minval = data[k]; }
      if ( data[k] > maxval ) { maxval = data[k]; }
    }

    /* The maxval is most probably set for unphysical regions so
       this in these cases the probabilities are set to zero  */
    for ( k = 0; k < totpix; k++ ) {
      if ( data[k] == maxval ) { data[k] = 0; }
      else { data[k] = exp( - (data[k] - minval) ); }
      totsum += data[k];
    }

    /* Normalize, so that the total probability equals 1 */
    for ( k = 0; k < totpix; k++ ) { data[k] /= totsum; }

    dex = (int *) malloc(totpix * sizeof(int));
    dexHeapSort( dex, data, totpix );

    /* Calculate the enclosed probabilities */
    totsum = 0;
    for ( k = totpix-1; k >= 0; k-- ) {
      totsum += data[dex[k]];
      data[dex[k]] = totsum;
    }

    /* Write the enclosed probability to a new fits extension */
    status = 0;
    fits_write_pix(outfptr, TFLOAT, firstpx , totpix, data, &status);
    if ( status != 0 ) {
      cerr << "Warning, FITSIO error: " << status << "\n";
    }
    fits_close_file(outfptr, &status);

    /*
    cout << "minval: " << minval << "\ntotsum: " << totsum << "\n";
    cout << "maxval: " << maxval << "\n";
    cout << "Number of axis: " << dim << naxes[0] << naxes[1] << "\n";
    */
  }

  exit(0);  
}


void dexHeapSort(int* dex, const float* data, int n)
{
  int i,j,k,top,cur;

  for (i=0 ; i<n ; ++i) dex[i]=i;

  if (n==1) return;
  
  k=(n/2);
  top=n-1;

  for (;;) {
    if (k) {                     // Still in heap creation phase
      --k;
      cur=dex[k];
    }
    else {                       // in heap-selection phase (l=0)
      cur=dex[top];              // clear a space at end of array
      dex[top]=dex[0];           // retire top of the heap into it
      if (--top == 0) {          // done with last promotion
        dex[0]=cur;              // lowest ranked value
        return;
      }
    }

    i=k;                         // cur has data[i] ; find its new spot
    j=k*2+1;                     // j points to first child of data[i]
    while (j<=top) {
      if (j<top && data[dex[j]]<data[dex[j+1]]) ++j; // pick the greatest child
      if (data[cur]<data[dex[j]]) {
        dex[i]=dex[j];                   // promote data[j] to slot i
        i=j;                               // then jump to j's children
        j=2*j+1;                           // to fill j's vacancy
      }
      else break;                // i is cur's level.
    }
    dex[i]=cur;                 // Put cur into its slot
  }
}

/* EOF */
