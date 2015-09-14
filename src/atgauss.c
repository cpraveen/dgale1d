#include "dg.h"
#include "dg1d.h"

void UatGauss(CELL * cell, REAL ** U)
{
   UINT iv, ig, ip;

   for(ig = 0; ig < cell->ng; ig++)
      for(iv = 0; iv < NVAR; iv++)
      {
         U[ig][iv] = 0.0;
         for(ip = 0; ip < cell->p; ip++)
            U[ig][iv] += cell->U[iv][ip] * ShapeFun(cell->xg[ig], cell, ip);
      }
}

// Compute solution at GLL nodes
void UatGLL(CELL * cell, REAL ** U)
{
   UINT iv, ig, ip;
   REAL xi, x;

   for(ig = 0; ig < cell->ngll; ig++)
   {
      xi = xgll[cell->ngll-1][ig];
      x  = 0.5*(1.0-xi)*cell->xl + 0.5*(1.0+xi)*cell->xr;
      for(iv = 0; iv < NVAR; iv++)
      {
         U[ig][iv] = 0.0;
         for(ip = 0; ip < cell->p; ip++)
            U[ig][iv] += cell->U[iv][ip] * ShapeFun(x, cell, ip);
      }
   }
}

// Compute solution at given value of x
void Uvect(CELL * cell, REAL x, REAL * U)
{
   UINT iv, ip;

   for(iv = 0; iv < NVAR; iv++)
   {
      U[iv] = 0.0;
      for(ip = 0; ip < cell->p; ip++)
         U[iv] += cell->U[iv][ip] * ShapeFun(x, cell, ip);
   }
}
