#include "dg.h"
#include "dg1d.h"

void Update(CELL * cell, FACE * face)
{
   UINT i, j, k;

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
      }
}
