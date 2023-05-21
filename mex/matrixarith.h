/*
 * Interface definition header file for matrixarith.c
 * Native C implementations of basic matrix arithmetic
 *
 * $Header: /home/mimo/CVSROOT/MIMO_SoC/SPHERE/matlab/Coded_OFDM_Simulation/matrixarith.h,v 1.2 2004/09/04 11:59:12 mimo Exp $
 *
 * Copyright (C) 2001 by Moritz Borgmann
 *
 * $Log: matrixarith.h,v $
 * Revision 1.2  2004/09/04 11:59:12  mimo
 * dos2unix
 *
 * Revision 1.1  2004/09/03 21:08:57  apburg
 * soft-viterbi
 *
 * Revision 1.1  2001/11/12 23:50:33  moriborg
 * Initial Checkin
 *
 *
 */

#ifndef __MATRIXARITH__
#define __MATRIXARITH__

#include "mexheader.h"

void matrixMultiply(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
void matrixSubtract(mxArrayConstPtr R, mxArrayConstPtr A, mxArrayConstPtr B);
double matrixFrobeniusNormSquared(mxArrayConstPtr A);

#endif
