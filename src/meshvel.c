#include <stdio.h>
#include <stdlib.h>
#include "dg.h"
#include "dg1d.h"

void MeshVel(FACE* face)
{
   UINT i;
   UINT siter=1; // no of smoothing iterations
   REAL UL[NVAR], UR[NVAR], UA[NVAR], vl, vr;
   
   if(ALE==0) return;
   
   face[0].w    = 0.0;
   face[NF-1].w = 0.0;
   
   for(i = 1; i < NF-1; i++)
      if(face[i].active)
      {
         Uvect(face[i].lcell, face[i].x, UL);
         Uvect(face[i].rcell, face[i].x, UR);
         
         vl = UL[1]/UL[0];
         vr = UR[1]/UR[0];
         face[i].w = 0.5 * (vl + vr);
//         RoeAverage(UL, UR, UA);
//         face[i].w = UA[1]/UA[0];
      }
   
   // smooth mesh velocity
   for(int it=0; it<siter; ++it)
   {
      for(i = 1; i < NF-1; i++)
         if(face[i].active)
         {
            face[i].w = (face[i].w + face[i].lcell->wl + face[i].rcell->wr)/3.0;
         }
      
      // copy to cells
      for(i = 1; i < NF-1; i++)
         if(face[i].active)
         {
            face[i].lcell->wr = face[i].w;
            face[i].rcell->wl = face[i].w;
         }
   }
}

// Do linear interpolation to compute mesh velocity at x
double GetMeshVel(CELL* cell, REAL x)
{
   return ((cell->xr - x) * cell->wl + (x - cell->xl) * cell->wr)/cell->h;
}