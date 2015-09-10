#ifndef _DG1D
#define _DG1D  1

#include <stdbool.h>

#define LINCON 1
#define BURGER 2
#define EULER  3

// Available test cases
#define SOD      0
#define BLAST    1
#define LAX      2
#define LOWD     3
#define SHUOSHER 4

// Type of bc
#define FREE      0
#define FIXED     1
#define PERIODIC  2

REAL mass0, mass1[2][2], mass2[3][3], mass3[4][4], mass4[5][5], mass5[6][6];
REAL cfl, dt, finaltime;
REAL XS;                        /* Shock position */
REAL xmin, xmax;
REAL d_left, u_left, p_left;
REAL d_right, u_right, p_right;
REAL Mfact;
UINT ALE;
REAL dxmin, dxmax;
UINT test_case;
UINT bc_left, bc_right;

struct CELL
{
   REAL x, xl, xr, wl, wr, h, *xg;
   UINT p, ng;
   REAL **U, **Re;
   struct CELL *lcell, *rcell;
   bool active;
};
typedef struct CELL CELL;

struct FACE
{
   REAL x; // location of face
   REAL w; // velocity of face
   CELL *lcell, *rcell;
   bool active;
};
typedef struct FACE FACE;

CELL* Init();
FACE* InitFaces(CELL*);
void TimeStep(CELL*);
void SaveSol(CELL*);
void Flux(CELL*, FACE*);
void Update(CELL*, FACE*);
void Project(CELL*);
void Result(CELL*);
void MeshVel(FACE*);

void GaussInit ();
void GaussPoints (CELL *);
REAL ShapeFun (REAL, CELL *, UINT);
REAL ShapeFunDeriv (REAL, CELL *, UINT);

void Uvect (CELL * cell, REAL x, REAL * U);
void get_predictor(CELL *cell, REAL x, REAL t, REAL * U);
void EulerFlux (REAL * U, REAL w, REAL * flux);
void Jacobian (REAL* U, REAL A[][3]);
void RoeFlux (REAL * Ul, REAL * Ur, REAL w, REAL * flux);
void LFFlux (REAL * Ul, REAL * Ur, REAL w, REAL * flux);
void ECUSPFlux (REAL * Ul, REAL * Ur, REAL * flux);
void HLLCFlux (REAL * Ul, REAL * Ur, REAL * flux);
void AUSMDVFlux (REAL * Ul, REAL * Ur, REAL * flux);
void LFCFlux (REAL * Ul, REAL * Ur, REAL * flux);

void EigMat (REAL *, REAL[][3], REAL[][3]);
void Multi (REAL[][3], REAL *);

void MoveGrid(CELL*, FACE*);
REAL GetMeshVel(CELL* cell, REAL x);

#endif
