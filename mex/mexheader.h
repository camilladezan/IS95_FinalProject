/*
 * C header file with several useful definitions for MATLAB MEX programming
 *
 * $Header: /home/mimo/CVSROOT/MIMO_SoC/SPHERE/matlab/Coded_OFDM_Simulation/mexheader.h,v 1.2 2004/09/04 11:59:11 mimo Exp $
 *
 * Copyright (C) 2001 by Moritz Borgmann
 *
 * $Log: mexheader.h,v $
 * Revision 1.2  2004/09/04 11:59:11  mimo
 * dos2unix
 *
 * Revision 1.1  2004/09/03 21:08:57  apburg
 * soft-viterbi
 *
 * Revision 1.2  2001/11/14 02:26:05  moriborg
 * Defined nil
 *
 * Revision 1.1  2001/11/12 23:50:35  moriborg
 * Initial Checkin
 *
 *
 */

#ifndef __MEXHEADER__
#define __MEXHEADER__

#include <mex.h>
#include <math.h>
#include <string.h>

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define NIL 0L
#define nil 0L

typedef mxArray * mxArrayPtr;
typedef const mxArray * mxArrayConstPtr;
typedef double * doublePtr;

#endif
