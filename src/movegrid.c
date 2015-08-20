#include <stdio.h>
#include <stdlib.h>
#include "dg.h"
#include "dg1d.h"

void MoveGrid(CELL* cell, FACE* face)
{
   UINT i;
   
   // Fixed boundaries, needs to be changed for periodic case
   
   for(i = 1; i < NF-1; i++)
      if(face[i].active)
      {
         face[i].x += face[i].w * dt;
  
         face[i].lcell->xr = face[i].x;
         face[i].rcell->xl = face[i].x;
      }
   
   dxmin =  1.0e20;
   dxmax = -1.0e20;
   for(i=0; i<NC; ++i)
      if(cell[i].active)
      {
         cell[i].h = cell[i].xr - cell[i].xl;
         cell[i].x = 0.5*(cell[i].xl + cell[i].xr);
         GaussPoints(&cell[i]);
         dxmin = MIN(dxmin, cell[i].h);
         dxmax = MAX(dxmax, cell[i].h);
      }
}
