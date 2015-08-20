#include <math.h>
#include "dg.h"
#include "dg1d.h"

// x is at time t_n
// t is relative to t_n
void get_predictor(CELL *cell, REAL x, REAL t, REAL * U)
{
   UINT i, j;
   REAL w, vx, Ux[NVAR], A[NVAR][NVAR];
   
   Uvect(cell, x, U);
   
   // compute Jacobian
   Jacobian (U, A);
   
   // mesh velocity at x
   w = GetMeshVel(cell, x);
   
   for(i=0; i<NVAR; ++i)
   {
      A[i][i] -= w;
      Ux[i] = 0.0;
   }
   
   // compute gradient of solution at x
   for(i = 0; i < cell->p; i++)
   {
      vx = ShapeFunDeriv(x, cell, i);
      for(j = 0; j < NVAR; j++)
         Ux[j] += vx * cell->U[j][i];
   }
   
   for(i=0; i<NVAR; ++i)
      for(j=0; j<NVAR; ++j)
         U[i] -= t * A[i][j] * Ux[j];
}