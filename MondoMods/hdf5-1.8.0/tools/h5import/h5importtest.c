/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
#include "H5private.h"

/*
 * Name:
 *      h5importtest
 *
 * Description:
 *      This program creates files  that can be
 *      used to test the h5import program.
 *
 */
/*
 * Define names for test files
 */

int
main(void)
{
    int         nrow = 3, ncol = 4, npln = 5;
    int         i, j, k;
    FILE       *sp;

    float     b32r3[5][3][4];
    float     row4[3], col4[4], pln4[5];
    float     rowo4 = (float)11.0e0, colo4 = (float)21.0e0, plno4 = (float)51.0e0;
    float     rowi4 = (float)1.0e0, coli4 = (float)2.0e0, plni4 = (float)5.0e0;

    int     b32i3[5][3][4];
    int     row4i[3], col4i[4], pln4i[5];
    int     rowo4i = (int)11 , colo4i = (int)21 , plno4i = (int)51 ;
    int     rowi4i = (int)1 , coli4i = (int)2 , plni4i = (int)5 ;

#ifndef _WIN32
    long_long     b64i2[3][4], b64i3[5][3][4];
    long_long     row4i64[3], col4i64[4], pln4i64[5];
    long_long     rowo4i64 = (long_long)11 , colo4i64 = (long_long)21 , plno4i64 = (long_long)51 ;
    long_long     rowi4i64 = (long_long)1 , coli4i64 = (long_long)2 , plni4i64 = (long_long)5 ;
#endif

    short     b16i3[5][3][4];
    short     row4i16[3], col4i16[4], pln4i16[5];
    short     rowo4i16 = (short)11 , colo4i16 = (short)21 , plno4i16 = (short)51 ;
    short     rowi4i16 = (short)1 , coli4i16 = (short)2 , plni4i16 = (short)5 ;

    char     b8i3[5][3][4];
    char     row4i8[3], col4i8[4], pln4i8[5];
    char     rowo4i8 = (char)11 , colo4i8 = (char)21 , plno4i8 = (char)51 ;
    char     rowi4i8 = (char)1 , coli4i8 = (char)2 , plni4i8 = (char)5 ;

    double     b64r3[5][3][4];
    double     row8[3], col8[4], pln8[5];
    double     rowo8 = 11.0e0, colo8 = 21.0e0, plno8 = 51.0e0;
	double     rowi8 = 1.0e0, coli8 = 2.0e0, plni8 = 5.0e0;


    /*
     * initialize the row, column, and plane vectors
     *
     * row values start at 11 and increment by 1 => 11, 12, 13
     * column values start at 21 and increment by 2 => 21, 23, 25, 27
     * plane values start at 51 and increment by 5 => 51, 56, 61, 66, 71
     */


    /*
     * build array elements - rank 2
     *
     * element value = sum of row value and col values
     */

    row4[0] = rowo4;
    col4[0] = colo4;
    pln4[0] = plno4;

    row8[0] = rowo8;
    col8[0] = colo8;
    pln8[0] = plno8;

    row4i[0] = rowo4i;
    col4i[0] = colo4i;
    pln4i[0] = plno4i;

#ifndef _WIN32
    row4i64[0] = rowo4i64;
    col4i64[0] = colo4i64;
    pln4i64[0] = plno4i64;
#endif

    row4i16[0] = rowo4i16;
    col4i16[0] = colo4i16;
    pln4i16[0] = plno4i16;

    row4i8[0] = rowo4i8;
    col4i8[0] = colo4i8;
    pln4i8[0] = plno4i8;

    for (i = 1; i < nrow; i++)
    {
        row4[i] = row4[i - 1] + rowi4;
        row8[i] = row8[i - 1] + rowi8;
        row4i[i] = row4i[i - 1] + rowi4i;
#ifndef _WIN32
	row4i64[i] = row4i64[i - 1] + rowi4i64;
#endif
	row4i16[i] = row4i16[i - 1] + rowi4i16;
	row4i8[i] = row4i8[i - 1] + rowi4i8;
    }

    for (j = 1; j < ncol; j++)
    {
          col4[j] = col4[j - 1] + coli4;
          col8[j] = col8[j - 1] + coli8;
	  col4i[j] = col4i[j - 1] + coli4i;
#ifndef _WIN32
          col4i64[j] = col4i64[j - 1] + coli4i64;
#endif
	  col4i16[j] = col4i16[j - 1] + coli4i16;
	  col4i8[j] = col4i8[j - 1] + coli4i8;
    }
    for (k = 1; k < npln; k++)
    {
          pln4[k] = pln4[k - 1] + plni4;
          pln8[k] = pln8[k - 1] + plni8;
	  pln4i[k] = pln4i[k - 1] + plni4i;
#ifndef _WIN32
    	  pln4i64[k] = pln4i64[k - 1] + plni4i64;
#endif
	  pln4i16[k] = pln4i16[k - 1] + plni4i16;
	  pln4i8[k] = pln4i8[k - 1] + plni4i8;
   }

   for (i = 0; i < nrow; i++)
   {
     for (j = 0; j < ncol; j++)
     {
#ifndef _WIN32
  	b64i2[i][j] = row4i64[i] + col4i64[j];
#endif
     }
   }

    /*
     * build array elements - rank 3
     *
     * element value = sum of row value, col, and plane values
     */

    for (i = 0; i < nrow; i++)
      {
          for (j = 0; j < ncol; j++)
            {
                for (k = 0; k < npln; k++)
                  {
                      b32r3[k][i][j] = row4[i] + col4[j] + pln4[k];
                      b64r3[k][i][j] = row8[i] + col8[j] + pln8[k];
	              b32i3[k][i][j] = row4i[i] + col4i[j] + pln4i[k];
#ifndef _WIN32
	     	      b64i3[k][i][j] = row4i64[i] + col4i64[j] + pln4i64[k];
#endif
		      b16i3[k][i][j] = row4i16[i] + col4i16[j] + pln4i16[k];
		      b8i3[k][i][j] = row4i8[i] + col4i8[j] + pln4i8[k];
                  }
            }
      }

    /*
     * binary 32-bit file - rank 2 & 3
     */

#ifndef UNICOS


    	sp = fopen("txtin16", "w");
	for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        {
          for (j = 0; j < ncol; j++)
              (void) fprintf(sp, "%10u", b16i3[k][i][j]);
          (void) fprintf(sp, "\n");
        }
    (void) fclose(sp);

      sp = fopen("txtin32", "w");
      for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        {
          for (j = 0; j < ncol; j++)
              (void) fprintf(sp, "%10d", b32i3[k][i][j]);
          (void) fprintf(sp, "\n");
        }
      (void) fclose(sp);

    sp = fopen("bin32", "w");
    for (k = 0; k < npln; k++)
    for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b32i3[k][i][j], sizeof(int), 1, sp);
    (void) fclose(sp);


    sp = fopen("buin32", "w");
    for (k = 0; k < npln; k++)
    for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b32i3[k][i][j], sizeof(unsigned int), 1, sp);
    (void) fclose(sp);

#ifndef _WIN32

    sp = fopen("bin64-2", "w");
    for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b64i2[i][j], sizeof(long_long), 1, sp);
    (void) fclose(sp);
#endif

    sp = fopen("bfp32", "w");
    for (k = 0; k < npln; k++)
        for (i = 0; i < nrow; i++)
            for (j = 0; j < ncol; j++)
                (void) fwrite((char *) &b32r3[k][i][j],
                              sizeof(float), 1, sp);
    (void) fclose(sp);

        sp = fopen("bin16", "w");
	for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b16i3[k][i][j], sizeof(short), 1,
                          sp);
    (void) fclose(sp);

        sp = fopen("buin16", "w");
	for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b16i3[k][i][j], sizeof(unsigned short), 1,
                          sp);
    (void) fclose(sp);

#ifndef _WIN32

        sp = fopen("bin64-3", "w");
	for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b64i3[k][i][j], sizeof(long_long), 1,
                          sp);
    (void) fclose(sp);
#endif

    sp = fopen("bin8", "w");
    for (k = 0; k < npln; k++)
    	for (i = 0; i < nrow; i++)
        for (j = 0; j < ncol; j++)
            (void) fwrite((char *) &b8i3[k][i][j], sizeof(char), 1,
                          sp);
    (void) fclose(sp);

#endif

    /*
     * binary 64-bit file - rank 2 & 3
     */

    sp = fopen("bfp64", "w");
    for (k = 0; k < npln; k++)
        for (i = 0; i < nrow; i++)
            for (j = 0; j < ncol; j++)
                (void) fwrite((char *) &b64r3[k][i][j],
                              sizeof(double), 1, sp);
    (void) fclose(sp);
    return (0);
}

