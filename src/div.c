#include <stdio.h>
#include <stdlib.h>
#include "dg.h"
#include "dg1d.h"

void Flux(CELL * cell, FACE * face)
{
   UINT i, j, k, l;
   REAL UL[NVAR], UR[NVAR], fl[NVAR], UG[NVAR], flg[NVAR], v, vx;

   for(i = 0; i < NC; i++)
      if(cell[i].active)
         for(j = 0; j < NVAR; j++)
            for(k = 0; k < cell[i].p; k++)
               cell[i].Re[j][k] = 0.0;

   /* Loop over cell faces and find flux, periodic bc */
   for(i = 0; i < NF; i++)
      if(face[i].active)
      {
         // mid-point time integration
         get_predictor(face[i].lcell, face[i].x, 0.5*dt, UL);
         get_predictor(face[i].rcell, face[i].x, 0.5*dt, UR);
         
         switch (FLUX)
         {
            case LF:
               LFFlux(UL, UR, face[i].w, fl);
               break;
            case ROE:
               RoeFlux(UL, UR, face[i].w, fl);
               break;
            default:
               printf("Error: Flux number %d not defined\n", FLUX);
               exit(0);
         }
         
         /* Add interface flux to the cells */
         
         if(i > 0)
            for(j = 0; j < NVAR; j++)
               for(k = 0; k < face[i].lcell->p; k++)
               {
                  v = ShapeFun(face[i].x, face[i].lcell, k);
                  face[i].lcell->Re[j][k] += fl[j] * v;
               }
         
         if(i < NF-1)
            for(j = 0; j < NVAR; j++)
               for(k = 0; k < face[i].rcell->p; k++)
               {
                  v = ShapeFun(face[i].x, face[i].rcell, k);
                  face[i].rcell->Re[j][k] -= fl[j] * v;
               }
      }

   /* Flux quadrature */
   for(i = 0; i < NC; i++)
      if(cell[i].active)
         for(j = 0; j < cell[i].ng; j++)
         {
            double w = GetMeshVel(&cell[i], cell[i].xg[j]);
            get_predictor(&cell[i], cell[i].xg[j], 0.5*dt, UG);
            EulerFlux(UG, w, flg);
            for(k = 0; k < cell[i].p; k++) {
               vx = ShapeFunDeriv(cell[i].xg[j], &cell[i], k);
               for(l = 0; l < NVAR; l++)
                  cell[i].Re[l][k] -=
                  0.5 * cell[i].h * flg[l] * vx * wg[cell[i].ng - 1][j];
            }
         }

}
