#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "dg.h"
#include "dg1d.h"

void Update(CELL * cell, FACE * face)
{
   UINT i, j, k;
   REAL p;

   for(i = 0; i < NC; i++)
      if(cell[i].active)
      {
         for(j = 0; j < NVAR; j++)
            for(k = 0; k < cell[i].p; k++)
               cell[i].U[j][k] = cell[i].h * cell[i].U[j][k] - dt * cell[i].Re[j][k];
      }
   
   MoveGrid(cell, face);
   
   for(i = 0; i < NC; i++)
      if(cell[i].active)
      {
         for(j = 0; j < NVAR; j++)
            for(k = 0; k < cell[i].p; k++)
               cell[i].U[j][k] /= cell[i].h;
         
         // Check positivity
         if(cell[i].U[0][0] < 0.0)
         {
            printf("Density is negative in cell %d !!!\n", i);
            exit(0);
         }
         p = (GAMMA-1.0)*(cell[i].U[2][0] - 0.5*pow(cell[i].U[1][0],2)/cell[i].U[0][0]);
         if(p < 0.0)
         {
            printf("Pressure is negative in cell %d !!!\n", i);
            exit(0);
         }
      }
}
