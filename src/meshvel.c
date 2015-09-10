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
   
   // First face
   i = 0;
   if(bc_left == FIXED)
   {
      face[i].w = 0.0;
   }
   else
   {
      Uvect(face[i].rcell, face[i].x, UR);
      face[i].w = UR[1]/UR[0];
   }
   face[i].rcell->wl = face[i].w;
   
   // Last face
   i = NF-1;
   if(bc_right == FIXED)
   {
      face[i].w = 0.0;
   }
   else
   {
      Uvect(face[i].lcell, face[i].x, UL);
      face[i].w = UL[1]/UL[0];
   }
   face[i].lcell->wr = face[i].w;
   
   // Interior faces
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
         
         face[i].lcell->wr = face[i].w;
         face[i].rcell->wl = face[i].w;
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