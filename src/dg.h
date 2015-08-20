#ifndef _DG

#define _DG         1

#define REAL        double
#define UINT        int

#define LF          1
#define ECUSP       2
#define HLLC        3
#define AUSMDV      4
#define LFC         5

#define CONSTANT    1
#define LINEAR      2
#define QUADRATIC   3
#define CUBIC       4
#define BIQUADRATIC 5

#define GAMMA       1.4
#define NVAR        3

#define MIN(a,b)   ( (a) < (b) ? (a) : (b) )
#define MAX(a,b)   ( (a) > (b) ? (a) : (b) )

/* NC = number of cells
 * NF = number of faces = NC + 1
 * NVAR = number of variables
 * NG = number of Gauss integration points
 */
UINT NC, NF, NG, PORD, FLUX, NPLT;

/* xg = Gauss integration points in [-1,+1]
 * wg = corresponding weights
 */
REAL xg[10][10], wg[10][10];

#endif
