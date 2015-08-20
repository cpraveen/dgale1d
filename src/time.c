#include <math.h>
#include "dg.h"
#include "dg1d.h"

void TimeStep(CELL * cell)
{
   UINT i;
   REAL d, u, p, c, dtc, w;
   REAL beta = 0.1;
   
   dt = 1.0e20;

   for(i = 0; i < NC; i++)
      if(cell[i].active)
      {
         d   = cell[i].U[0][0];
         u   = cell[i].U[1][0] / d;
         p   = (GAMMA - 1.0) * (cell[i].U[2][0] - 0.5 * d * u * u);
         c   = sqrt(GAMMA * p / d);
         w   = 0.5 * (cell[i].wl + cell[i].wr);
         dtc = cfl * cell[i].h / (fabs(u-w) + c);
         dt  = MIN(dtc, dt);
         
         dtc = beta * cell[i].h / (fabs(cell[i].wr-cell[i].wl) + 1.0e-14);
         dt  = MIN(dtc, dt);
      }

}
